#!/usr/bin/env python3

from PIL import Image
from glob import glob
import cv2
import os

for i in range(6):
    if not os.path.isdir("cand_grade_{0}".format(i)):
        os.mkdir("cand_grade_{0}".format(i))


images = glob('*png')
print(images[0])
for i in images:
    """
    im = Image.open(i)
    im.show()
    cand_grade = input()
    print(cand_grade)
    im.close()
    """
    print(i)
    im = cv2.imread(i)
    cv2.namedWindow('image',cv2.WINDOW_NORMAL)
    cv2.resizeWindow('image', 1400, 1800)
    cv2.imshow('image', im)
    cv2.waitKey(3000)
    cv2.destroyAllWindows()
    cand_grade = input()
    os.rename(i, "cand_grade_{0}/{1}".format(cand_grade, i))

