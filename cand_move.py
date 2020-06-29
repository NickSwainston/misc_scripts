#!/usr/bin/env python3

try:
    from PIL import Image
except ImportError:
    import Image
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
    cv2.waitKey(2000)
    cv2.destroyAllWindows()
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

