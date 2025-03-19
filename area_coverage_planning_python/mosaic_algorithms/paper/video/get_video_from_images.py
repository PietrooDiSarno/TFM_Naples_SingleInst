import cv2
import os
import re


def create_video_from_images(image_folder, fps=30, target_resolution=(1920, 1080)):
    # Ensure the output video is saved in the same folder as the images
    output_video_path = os.path.join(image_folder, "output_video.mp4")

    # Get the list of images in the folder
    images = [f for f in os.listdir(image_folder) if f.lower().endswith(('png', 'jpg', 'jpeg'))]

    # Sort images based on the numerical value in the filenames
    def extract_number(filename):
        match = re.search(r'(\d+)', filename)
        return int(match.group(1)) if match else float('inf')  # Use 'inf' for filenames without numbers

    images.sort(key=extract_number)

    if not images:
        print("No images found in the specified folder.")
        return

    # Define the codec and create VideoWriter object
    fourcc = cv2.VideoWriter_fourcc(*'mp4v')
    video = cv2.VideoWriter(output_video_path, fourcc, fps, target_resolution)

    # Add each image to the video
    for image_name in images:
        image_path = os.path.join(image_folder, image_name)
        frame = cv2.imread(image_path)
        if frame is None:
            print(f"Warning: Could not read {image_path}, skipping.")
            continue

        # Resize frame to target resolution
        frame_resized = cv2.resize(frame, target_resolution)
        video.write(frame_resized)

    # Release the video writer
    video.release()
    print(f"Video saved to {output_video_path}")


# Usage
image_folder = "C:/Users/pietr/PycharmProjects/TFM_Naples_SingleInst/area_coverage_planning_python/mosaic_algorithms/paper/video/png_images"  # Replace with the path to your folder containing images
create_video_from_images(image_folder)