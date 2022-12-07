#include <opencv2/opencv.hpp>
#include <opencv2/highgui.hpp>
#include <iostream>
#include <vector>

// g++ `pkg-config --cflags opencv4` cv2.cpp `pkg-config --libs opencv4`

using namespace std;
using namespace cv;

int main(int argc, char **argv)
{
    std::vector<int> lower = {100,100,125};
    std::vector<int> upper = {255,255,255};
    VideoCapture cap1(0);
    if (!cap1.isOpened())
        return -1;

    namedWindow("pierwsze", WINDOW_AUTOSIZE);
    namedWindow("detected", WINDOW_AUTOSIZE);

    while (true)
    {
        Mat f1, dst, detected,dilated;
        cap1.read(f1);
        cvtColor(f1, dst, COLOR_BGR2HSV);

        inRange(dst, lower,upper, detected);
        auto kernel = getStructuringElement(MORPH_ELLIPSE,Size{5,5});
        erode(detected, dilated, kernel);
        dilate(dilated, dilated, kernel);

        for (int y = 0; y < dilated.rows; ++y){
            for (int x = 0; x < dilated.cols; ++x) {
                if(dilated.at<Vec3b>(y,x) == [0,0,0])
        }

        /*for (int y = 0; y < input_bgra.rows; ++y) {
            cv::Vec4b &pixel = input_bgra.at<cv::Vec4b>(y, x)
            for (int x = 0; x < input_bgra.cols; ++x) {
                cv::Vec4b &pixel = input_bgra.at<cv::Vec4b>(y, x);

                if (pixel[0] = 1) {
                    pixel[1] = 122;
                    pixel[2] = 122;
                    pixel[3] = 122;
                }
            }
        }
         */

        imshow("dilated", dilated);
        if (waitKey(1) == 27)
            break;


    }

    return 0;
}