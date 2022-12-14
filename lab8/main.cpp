#include <iostream>
#include <opencv2/highgui.hpp>
#include <opencv2/opencv.hpp>
#include <vector>

using namespace std;
using namespace cv;

bool compareContours(vector<Point> c1, vector<Point> c2) {
    double i = fabs(contourArea(Mat(c1)));
    double j = fabs(contourArea(Mat(c2)));
    return (i < j);
}
std::vector<int> lower = {100,100,100};
std::vector<int> upper = {255, 255, 255};

int main(int argc, char **argv) {

    VideoCapture camera(0);
    if (!camera.isOpened())
        return -1;

    while (waitKey(1) != 27) {
        Mat img;
        camera.read(img);

        Mat base_img;
        Mat blue_img;

        cvtColor(img, base_img, COLOR_BGR2HSV);
        inRange(base_img, Scalar(lower[0], lower[1], lower[2]), Scalar(upper[0], upper[1], upper[2]), blue_img);

        Mat kernel = getStructuringElement(MORPH_ELLIPSE, {1, 1});
        morphologyEx(blue_img, blue_img, MORPH_CLOSE, kernel);

        Canny(blue_img, blue_img, 100, 100);

        auto mat_kernel = getStructuringElement(MORPH_ELLIPSE, {50, 50});
        morphologyEx(blue_img, blue_img, MORPH_CLOSE, mat_kernel);

        vector<vector<Point>> conts;
        findContours(blue_img, conts, RETR_LIST, CHAIN_APPROX_SIMPLE);
        //drawContours(img, conts, -1, {0, 255, 0}, 2);

        sort(conts.begin(), conts.end(), compareContours);

        if (conts.size() > 1) {
            vector<Point> first = conts[conts.size() - 1];
            vector<Point> second = conts[conts.size() - 2];

            auto momentsFirst = moments(first, false);
            auto momentsSecond = moments(second, false);

            Point p = {(int) (momentsFirst.m10 / momentsFirst.m00), (int) (momentsFirst.m01 / momentsFirst.m00)};
            Point p_2 = {(int) ((momentsFirst.m10 / momentsFirst.m00) - 5), (int) ((momentsFirst.m01 / momentsFirst.m00) - 5)};
            Point p_3 = {(int) ((momentsFirst.m10 / momentsFirst.m00) - 10), (int) ((momentsFirst.m01 / momentsFirst.m00) - 10)};

            Point p1 = {(int) (momentsSecond.m10 / momentsSecond.m00), (int) (momentsSecond.m01 / momentsSecond.m00)};
            Point p1_2 = {(int) ((momentsSecond.m10 / momentsSecond.m00) - 5), (int) ((momentsSecond.m01 / momentsSecond.m00) - 5)};
            Point p1_3 = {(int) ((momentsSecond.m10 / momentsSecond.m00) - 10), (int) ((momentsSecond.m01 / momentsSecond.m00) - 10)};


            line(img, p, p1, cv::Scalar(0, 255, 0), 3, cv::LINE_8);
            line(img, p_2, p1_2, cv::Scalar(255, 255, 0), 3, cv::LINE_8);
            line(img, p_3, p1_3, cv::Scalar(0, 255, 255), 3, cv::LINE_8);

        }

        imshow("Draw 3 Lines", img);
    }
    return 0;
}