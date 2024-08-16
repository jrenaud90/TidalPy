
#include "material_.hpp"


void MaterialCC::build_shear_rheology(int rheology_model_int, const double* rheology_parameters)
{
    this->shear_rheology_model_int = rheology_model_int;
    switch (this->shear_rheology_model_int)
    {
    case ELASTIC_RHEOLOGY_INT:
        this->shear_rheology_ptr = std::make_shared<Elastic>(rheology_parameters);
        break;
    case NEWTON_RHEOLOGY_INT:
        this->shear_rheology_ptr = std::make_shared<Newton>(rheology_parameters);
        break;
    case MAXWELL_RHEOLOGY_INT:
        this->shear_rheology_ptr = std::make_shared<Maxwell>(rheology_parameters);
        break;
    case VOIGT_RHEOLOGY_INT:
        this->shear_rheology_ptr = std::make_shared<Voigt>(rheology_parameters);
        break;
    case BURGERS_RHEOLOGY_INT:
        this->shear_rheology_ptr = std::make_shared<Burgers>(rheology_parameters);
        break;
    case ANDRADE_RHEOLOGY_INT:
        this->shear_rheology_ptr = std::make_shared<Andrade>(rheology_parameters);
        break;
    case SUNDBERG_RHEOLOGY_INT:
        this->shear_rheology_ptr = std::make_shared<SundbergCooper>(rheology_parameters);
        break;
    default:
        printf("Unsupported rheology model int %d encountered in Material class.", rheology_model_int);
        break;
    }
}

void MaterialCC::build_bulk_rheology(int rheology_model_int, const double* rheology_parameters)
{
    this->bulk_rheology_model_int = rheology_model_int;
    switch (this->bulk_rheology_model_int)
    {
    case ELASTIC_RHEOLOGY_INT:
        this->bulk_rheology_ptr = std::make_shared<Elastic>(rheology_parameters);
        break;
    case NEWTON_RHEOLOGY_INT:
        this->bulk_rheology_ptr = std::make_shared<Newton>(rheology_parameters);
        break;
    case MAXWELL_RHEOLOGY_INT:
        this->bulk_rheology_ptr = std::make_shared<Maxwell>(rheology_parameters);
        break;
    case VOIGT_RHEOLOGY_INT:
        this->bulk_rheology_ptr = std::make_shared<Voigt>(rheology_parameters);
        break;
    case BURGERS_RHEOLOGY_INT:
        this->bulk_rheology_ptr = std::make_shared<Burgers>(rheology_parameters);
        break;
    case ANDRADE_RHEOLOGY_INT:
        this->bulk_rheology_ptr = std::make_shared<Andrade>(rheology_parameters);
        break;
    case SUNDBERG_RHEOLOGY_INT:
        this->bulk_rheology_ptr = std::make_shared<SundbergCooper>(rheology_parameters);
        break;
    default:
        printf("Unsupported rheology model int %d encountered in Material class.", rheology_model_int);
        break;
    }
}