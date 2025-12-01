from PIL import Image, ImageChops

def trim(img):
    # Detect background color
    bg = Image.new(img.mode, img.size, img.getpixel((0,0)))
    
    # Find differences between the image and the background
    diff = ImageChops.difference(img, bg)
    
    # Convert difference to greyscale to compress channels
    diff = diff.convert("L")
    
    # Create bounding box of non-background content
    bbox = diff.getbbox()
    
    # Only crop if there is something to crop
    return img.crop(bbox) if bbox else img

def stack(*img_names: tuple[str], top_bottom=True):

    buffer = 0.05 # fraction of image size to use as buffer

    images = list()
    max_width = 0
    total_width = 0
    max_height = 0
    total_height = 0
    for img_name in img_names:
        images.append(trim(Image.open(img_name)))
        max_width = max(max_width, images[-1].width)
        max_height = max(max_height, images[-1].height)
    for img in images:
        if top_bottom:
            img.resize((max_width, img.height))
        else:
            img.resize((img.width, max_height))
        total_width += images[-1].width
        total_height += images[-1].height

    if top_bottom:
        # Create a new blank image with combined height
        buffer_size = int(max_height * buffer)
        total_buffer = buffer_size * (len(images) - 1)
        combined = Image.new("RGB", (max_width, total_height + total_buffer))

        current_y = 0
        for i, img in enumerate(images):
            combined.paste(img, (0, current_y))
            current_y += img.height + buffer_size
    else:
        # Create a new blank image with combined height
        buffer_size = int(max_width * buffer)
        total_buffer = buffer_size * (len(images) - 1)
        combined = Image.new("RGB", (total_width + total_buffer, max_height))

        current_x = 0
        for i, img in enumerate(images):
            combined.paste(img, (current_x, 0))
            current_x += img.width + buffer_size

    return combined


combos = [
    ("trappist1e_2d_heating_tidally_locked.png", "trappist1e_2d_heating_3-2_spin-orbit_resonance.png", "trappist1e_2d_heating_3-2_sor_15deg_obliquity.png"),
]

for i, img_names in enumerate(combos):
    combined_image = stack(*img_names, top_bottom=True)
    output_name = f"stacked_images_set_{i}.png"
    combined_image.save(output_name)
    # combined_image.show()
