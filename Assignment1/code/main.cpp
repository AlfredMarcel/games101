#include "Triangle.hpp"
#include "rasterizer.hpp"
#include <eigen3/Eigen/Eigen>
#include <iostream>
#include <opencv2/opencv.hpp>
#include <cmath>

constexpr double MY_PI = 3.1415926;

Eigen::Matrix4f get_view_matrix(Eigen::Vector3f eye_pos)
{
    Eigen::Matrix4f view = Eigen::Matrix4f::Identity();

    Eigen::Matrix4f translate;
    translate << 1, 0, 0, -eye_pos[0], 0, 1, 0, -eye_pos[1], 0, 0, 1,
        -eye_pos[2], 0, 0, 0, 1;

    view = translate * view;

    return view;
}

Eigen::Matrix4f get_model_matrix(float rotation_angle)
{
    Eigen::Matrix4f model = Eigen::Matrix4f::Identity();

    float rot = rotation_angle/180.0*MY_PI;

    Eigen::Matrix4f translate;
    // z
    translate << cos(rot), -sin(rot), 0, 0, sin(rot), cos(rot), 0, 0,
         0, 0, 1, 0, 0, 0, 0, 1; 

    // x
    //translate << 1,0,0,0,0,cos(rot),-sin(rot),0,0,sin(rot),cos(rot),0,0,0,0,1;

    // y 
    //translate << cos(rot),0,sin(rot),0,0,1,0,0,-sin(rot),0,cos(rot),0,0,0,0,1;

    model = translate * model;

    // TODO: Implement this function
    // Create the model matrix for rotating the triangle around the Z axis.
    // Then return it.
    
    return model;
}

Eigen::Matrix4f get_rotation(Vector3f axis, float angle){
    float rot = angle/180*MY_PI;
    Eigen::Matrix3f tmp;
    tmp<< 0,-axis[2],axis[1],axis[2],0,-axis[0],-axis[1],axis[0],0;

    Eigen::Matrix3f rotation = Eigen::Matrix3f::Identity();
    Eigen::Matrix4f res = Eigen::Matrix4f::Identity();

    rotation = cos(rot)*rotation+(1-cos(rot))*axis*axis.transpose()+sin(rot)*tmp;
    res.block<3,3>(0,0) = rotation;

    return res;
}

Eigen::Matrix4f get_projection_matrix(float eye_fov, float aspect_ratio,
                                      float zNear, float zFar)
{
    // Students will implement this function

    Eigen::Matrix4f projection = Eigen::Matrix4f::Identity();
    float fy = 1.0/tan(eye_fov/2/180*MY_PI);
    float fx = fy/aspect_ratio;
    float c = -(zNear+zFar)/(zFar-zNear);
    float d = -2*zNear*zFar/(zFar-zNear);

    Eigen::Matrix4f translate;
    translate << fx, 0, 0, 0, 0, fy, 0, 0, 0, 0, c, d, 0, 0, -1, 0; 

    projection = translate * projection;
    return projection;

    // TODO: Implement this function
    // Create the projection matrix for the given parameters.
    // Then return it.

}

int main(int argc, const char** argv)
{
    float angle = 0;
    bool command_line = false;
    std::string filename = "output.png";

    if (argc >= 3) {
        command_line = true;
        angle = std::stof(argv[2]); // -r by default
        if (argc == 4) {
            filename = std::string(argv[3]);
        }
        else
            return 0;
    }

    rst::rasterizer r(700, 700);

    Eigen::Vector3f eye_pos = {0, 0, 5};

    std::vector<Eigen::Vector3f> pos{{2, 0, -2}, {0, 2, -2}, {-2, 0, -2}};

    std::vector<Eigen::Vector3i> ind{{0, 1, 2}};

    Eigen::Vector3f axis = {1,1,1};

    auto pos_id = r.load_positions(pos);
    auto ind_id = r.load_indices(ind);

    int key = 0;
    int frame_count = 0;

    if (command_line) {
        r.clear(rst::Buffers::Color | rst::Buffers::Depth);

        r.set_model(get_model_matrix(angle));
        r.set_view(get_view_matrix(eye_pos));
        r.set_projection(get_projection_matrix(45, 1, 0.1, 50));

        r.draw(pos_id, ind_id, rst::Primitive::Triangle);
        cv::Mat image(700, 700, CV_32FC3, r.frame_buffer().data());
        image.convertTo(image, CV_8UC3, 1.0f);

        cv::imwrite(filename, image);

        return 0;
    }

    while (key != 27) {
        r.clear(rst::Buffers::Color | rst::Buffers::Depth);

        // r.set_model(get_model_matrix(angle));
        r.set_model(get_rotation(axis,angle));
        r.set_view(get_view_matrix(eye_pos));
        r.set_projection(get_projection_matrix(45, 1, 0.1, 50));

        r.draw(pos_id, ind_id, rst::Primitive::Triangle);

        cv::Mat image(700, 700, CV_32FC3, r.frame_buffer().data());
        image.convertTo(image, CV_8UC3, 1.0f);
        cv::imshow("image", image);
        key = cv::waitKey(10);

        std::cout << "frame count: " << frame_count++ << '\n';

        if (key == 'a') {
            angle += 10;
        }
        else if (key == 'd') {
            angle -= 10;
        }

        if(key == 'w'){
            eye_pos[2] +=1;
        }else if(key == 's'){
            eye_pos[2] -=1;
        }
    }

    return 0;
}
