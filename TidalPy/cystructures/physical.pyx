# distutils: language = c
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False

from libc.math cimport NAN, isnan, M_PI, fabs

# Constants
cdef double PI_4 = M_PI * 4.
cdef double PI_4_3 = PI_4 / 3.

# Epsilons for comparisons
cdef double EPS_RADIAL = 1.0e-6  # Any radial differences smaller than this will be ignored.


cdef class PhysicalStructure:
    
    def __init__(
        self,
        str name,
        double radius_upper = NAN,
        double radius_lower = NAN,
        double thickness = NAN,
        bint auto_update_physical = True
        ):

        # Set basic information
        self.name = name
        self.error_code = 0

        # Put placeholders on other properties
        self._radius_upper = NAN
        self._radius_lower = NAN
        self._thickness = NAN
        self._volume = NAN
        self._surface_area_lower = NAN
        self._surface_area_upper = NAN

        # Update physical properties
        if auto_update_physical:
            self.error_code = self.cf_change_radials(
                new_radius_lower = radius_lower,
                new_radius_upper = radius_upper,
                new_thickness = thickness,
                auto_update_volume_area=True
                )
        if self.error_code == -1:
            # -1 == ValueError
            raise ValueError(f"Error found in {self.name}:\n\t{self.error_message}.")
    
    cdef void _log(
            self,
            char* message,
            bint error_code = 0
            ) noexcept nogil:
        # TODO: Replace with actual logger.
        if error_code >= 0:
            self._status_message = message
        elif error_code == -1:
            self._error_message = message
            self.error_code = error_code
    
    cdef int cf_change_radials(
            self,
            double new_radius_lower = NAN,
            double new_radius_upper = NAN,
            double new_thickness = NAN,
            bint auto_update_volume_area = True
            ) noexcept nogil:
        """Change one or more of the radial properties."""

        cdef int error_code = 0
        cdef bint something_changed = False
        cdef bint radius_upper_provided = False
        cdef bint radius_lower_provided = False
        cdef bint thickness_provided = False
        
        # Parse input check for basic value errors.
        if not isnan(new_radius_upper):
            radius_upper_provided = True
            if new_radius_upper < 0.:
                self._log("(cf_change_radials):: negative value provided for `radius_upper`.", -1)
                return -1
        
        if not isnan(new_radius_lower):
            radius_lower_provided = True
            if new_radius_lower < 0.:
                self._log("(cf_change_radials):: negative value provided for `radius_lower`.", -1)
                return -1
        
        if not isnan(new_thickness):
            thickness_provided = True
            if new_thickness < 0.:
                self._log("(cf_change_radials):: negative value provided for `thickness`.", -1)
                return -1

        if thickness_provided and radius_lower_provided and radius_upper_provided:
            # All three provided. Make sure they are consistent.
            if new_thickness != (new_radius_upper - new_radius_lower):
                self._log("(cf_change_radials):: `radius_lower`, `radius_upper`, and `thickness` provided; however values found to be inconsistent with each other.", -2)
                return -2
        
        # Prioritize on the lower radius as this is least likely to change.
        if radius_lower_provided:
            if isnan(self._radius_lower):
                self._radius_lower = new_radius_lower
                something_changed = True
            else:
                # Radius lower already set. See if the new value is significantly different.
                if fabs(self._radius_lower - new_radius_lower) > EPS_RADIAL:
                    self._radius_lower = new_radius_lower
                    something_changed = True
        
        if radius_upper_provided:
            if isnan(self._radius_upper):
                self._radius_upper = new_radius_upper
                something_changed = True
            else:
                # Radius upper already set. See if the new value is significantly different.
                if fabs(self._radius_upper - new_radius_upper) > EPS_RADIAL:
                    self._radius_upper = new_radius_upper
                    something_changed = True
            
        if thickness_provided:
            if isnan(self._thickness):
                self._thickness = new_thickness
                something_changed = True
            else:
                # Thickness already set. See if the new value is significantly different.
                if fabs(self._thickness - new_thickness) > EPS_RADIAL:
                    self._thickness = new_thickness
                    something_changed = True
            
        # Change other radials based on these changes.
        if something_changed:
            if radius_lower_provided:
                if thickness_provided:
                    self._radius_upper = self._radius_lower + self._thickness
                elif radius_upper_provided:
                    self._thickness = self._radius_upper - self._radius_lower
                else:
                    # Only radius lower changed. If thickness is provided; update radius upper; otherwise update thickness.
                    if not isnan(self._thickness):
                        self._radius_upper = self._radius_lower + self._thickness
                    elif not isnan(self._radius_upper):
                        self._thickness = self._radius_upper - self._radius_lower
            if radius_upper_provided:
                if thickness_provided:
                    self._radius_lower = self._radius_upper - self._thickness
                elif radius_lower_provided:
                    # Covered by previous case.
                    pass
                else:
                    # Only radius upper changed. Assume that radius lower has not changed. Change thickness first if possible.
                    if not isnan(self._radius_lower):
                        self._thickness = self._radius_upper - self._radius_lower
                    elif not isnan(self._thickness):
                        self._radius_lower = self._radius_upper - self._thickness
            if thickness_provided:
                # Provided cases are covered with the above. Skip right to the presets.
                if not isnan(self._radius_lower):
                    self._radius_upper = self._radius_lower + self._thickness
                elif not isnan(self._radius_upper):
                    self._radius_lower = self._radius_upper - self._thickness
            
            # Check for bad values.
            if self._thickness < 0.:
                self._log("(cf_change_radials):: Negative `thickness` encountered after updating with `radius_upper` and/or `radius_lower`.", -1)
                return -1
            if self._radius_lower > self._radius_upper:
                self._log("(cf_change_radials):: Lower radius found to be larger than upper radius.", -1)
                return -1
            elif self._radius_lower < 0.:
                self._log("(cf_change_radials):: Negative lower radius encountered.", -1)
                return -1
            if self._radius_upper < 0.:
                self._log("(cf_change_radials):: Negative upper radius encountered.", -1)
                return -1

            # Call update functions 
            if auto_update_volume_area:
                error_code = self.cf_update_volume_area()
                if error_code != 0:
                    return error_code

        else:
            self._log("(cf_change_radials):: Called but no changed were made.", 2)

        return 0

    cdef int cf_update_volume_area(
            self
            ) noexcept nogil:
        """Update volumes and areas using state properties."""
        
        cdef double r2_upper = self._radius_upper**2
        cdef double r2_lower = self._radius_lower**2

        self._volume = PI_4_3 * ((self._radius_upper * r2_upper) - (self._radius_lower * r2_lower))
        self._surface_area_upper = PI_4 * r2_upper
        self._surface_area_lower = PI_4 * r2_lower

        return 0
    
    def change_radials(
            self,
            double new_radius_upper = NAN,
            double new_radius_lower = NAN,
            double new_thickness = NAN,
            bint auto_update_volume_area = True
            ):
        
        self.error_code = \
            self.cf_change_radials(
                new_radius_upper, new_radius_lower, new_thickness,
                auto_update_volume_area
                )
        
        if self.error_code == -1:
            raise ValueError(f"Error in {self.name}:\n\t{self.error_message}.")
        elif self.error_code == -2:
            raise AttributeError(f"Error in {self.name}:\n\t{self.error_message}.")


    # Cython defined properties
    @property
    def radius_upper(self):
        return self._radius_upper

    @property
    def radius_lower(self):
        return self._radius_lower
    
    @property
    def thickness(self):
        return self._thickness

    @property
    def volume(self):
        return self._volume
    
    @property
    def surface_area_lower(self):
        return self._surface_area_lower

    @property
    def surface_area_upper(self):
        return self._surface_area_upper

    # Aliased Properties
    @property
    def radius(self):
        return self._radius_upper
    
    @property
    def surface_area(self):
        return self._surface_area_upper
    
    # Derived Properties    
    @property
    def error_message(self):
        return str(self._error_message, encoding="UTF-8")
    
    @property
    def status_message(self):
        return str(self._status_message, encoding="UTF-8")
