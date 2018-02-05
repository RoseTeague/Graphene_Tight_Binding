"""Module for solver of with TB split operator
"""
import numpy as np
import matplotlib as mpl
import cv2
import os


def MakeMovie(title = 'Tight Binding Video'):

    image_folder = 'Images'
    video_name = title+'.mp4'

    #images = [img for img in os.listdir(image_folder) if img.endswith(".png")]


    images = []
    for f in os.listdir(image_folder):
        if f.endswith('.png'):
            images.append(int(f.strip('.png')))
    images.sort()

    frame = cv2.imread(os.path.join(image_folder, str(images[0])+'.png'))
    height, width, layers = frame.shape

    video = cv2.VideoWriter(video_name, -1, 10, (width,height),)

    for image in images:
        video.write(cv2.imread(os.path.join(image_folder, str(image)+'.png')))
        os.remove('Images/'+str(image)+'.png')

    cv2.destroyAllWindows()
    video.release()
