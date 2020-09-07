#!/usr/bin/env python3

import argparse
from glob import glob
import os

parser = argparse.ArgumentParser(description="""Simple graphical interface for candidate sorting""")
parser.add_argument("-m", "--method", type=str, default='cv2',
        help='The method to display the image. Either "matplotlib" or "cv2". Default: cv2.')
args = parser.parse_args()

if args.method == "cv2":
    try:
        from PIL import Image
    except ImportError:
        import Image
    import cv2
elif args.method == "matplotlib":
    import matplotlib.pyplot as plt
    import matplotlib.image as mpimg
    fig=plt.figure(figsize=(1250, 1000), dpi=1)
    plt.ion()
    plt.axis('off')
    fig.set_tight_layout(True)
    plt.show()
else:
    print("Incorrect method given")
    exit()

images = glob('*png')
for i in range(6):
    if not os.path.isdir("cand_grade_{0}".format(i)):
        os.mkdir("cand_grade_{0}".format(i))

for i in images:
    print(i)
    if args.method == "cv2":
        im = cv2.imread(i)
        cv2.namedWindow('image',cv2.WINDOW_NORMAL)
        cv2.resizeWindow('image', 1400, 1800)
        cv2.imshow('image', im)
        cv2.waitKey(2000)
        cv2.destroyAllWindows()
    elif args.method == "matplotlib":
        image=mpimg.imread(i)
        imgplot = plt.imshow(image)
        plt.show()
    cand_grade = input()
    os.rename(i, "cand_grade_{0}/{1}".format(cand_grade, i))
    pfd = i[:-4]
    tar = pfd + ".tar.gz"
    if os.path.isfile(pfd):
        os.rename(pfd, "cand_grade_{0}/{1}".format(cand_grade, pfd))
    elif "singlepulse" not in pfd:
        print('WARNING no pfd file')
    if os.path.isfile(tar):
        os.rename(tar, "cand_grade_{0}/{1}".format(cand_grade, tar))
    bestprof = pfd + '.bestprof'
    if os.path.isfile(bestprof):
        os.rename(bestprof, "cand_grade_{0}/{1}".format(cand_grade, bestprof))

