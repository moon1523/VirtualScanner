/*
 * Software License Agreement (BSD License)
 *
 *  Copyright (c) 2010, Willow Garage, Inc.
 *  All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions
 *  are met:
 *
 *   * Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *   * Redistributions in binary form must reproduce the above
 *     copyright notice, this list of conditions and the following
 *     disclaimer in the documentation and/or other materials provided
 *     with the distribution.
 *   * Neither the name of the copyright holder(s) nor the names of its
 *     contributors may be used to endorse or promote products derived
 *     from this software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 *  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 *  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *  COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 *  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 *  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 *  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 *  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 *  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 *  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 *
 * $Id$
 *
 */

/**
  * \author Radu Bogdan Rusu
  *
  * @b virtual_scanner takes in a .ply or a .vtk file of an object model, and virtually scans it
  * in a raytracing fashion, saving the end results as PCD (Point Cloud Data) files. In addition,
  * it noisifies the PCD models, and downsamples them.
  * The viewpoint can be set to 1 or multiple views on a sphere.
  */
#include <string>
#include "boost.h"
#include "vtkGeneralTransform.h"
#include "vtkTransform.h"
#include "vtkCellLocator.h"

#include "functions.cpp"

using namespace pcl;

#define EPS 0.00001

const double pi = 3.14159265358979f;
struct ScanParameters
{
    double fov_vert;          // vertical fov in radian
    double fov_hor;           // horizontal fov in radian
    int    res_vert;          // vertical resolution
    int    res_hor;           // horizontal  resolution
    double distance;          // distance in mm
    double screen_dist;       // screen distance in mm
    int    rotInterval;       // rotation interval in degree
};
using namespace std;
void PrintUsage(){
    cout<<"Usage: ./scanner [options] model.ply viewpoints.ply"<<endl;
    cout<<"[options]"<<endl;
    cout<<"\t-dist <double> : scan distance in mm (default:3000)"<<endl;
    cout<<"                   or put three numbers for multiple distances [start] [interval] [#]"<<endl;
    cout<<"\t-rot  <int>    : rotation interval in degree (default:10)"<<endl;
    cout<<"\t-path <string> : directory path for result files (default: .)"<<endl;
    cout<<"\t-cv   [depth/color/both] : print openCV matrix in YML"<<endl;
    cout<<"\t-ply  [0/1]    : print point cloud w/ color in PLY [0:object fixed/1:camera fixed]"<<endl;
    cout<<"\t-pcd  [0/1]    : quick print of point cloud w/o color in PCD [0:object fixed/1:camera fixed]"<<endl;
    exit(1);
}
void PrintKeyUsage(){
    cout<<"[Key Usage]"<<endl;
    cout<<"\tw : write template"<<endl;
    cout<<"\tv : to next view"<<endl;
    cout<<"\tf : show the number of failed templates"<<endl;
}

void EstimateColor(vtkPolyData* data, vtkIdType cellId, double* point, double* rgb);
void drawResponse(const std::vector<cv::linemod::Template>& templates,
    int num_modalities, cv::Mat& dst, cv::Point offset, int T);
static void writeLinemod(const cv::Ptr<cv::linemod::Detector>& detector, const std::string& filename);

int main (int argc, char** argv)
{
    //arguements (file names)
    if(argc<3) PrintUsage();
    string modelFileName = argv[argc-2];
    string viewFileName  = argv[argc-1];

    //set scan parameters
    ScanParameters scan;
    scan.fov_hor  = 90.*pi/180.;
    scan.fov_vert = 59.*pi/180.;
    scan.res_hor  = 1280;
    scan.res_vert  = 720;
    scan.distance = 3000; //mm
    scan.rotInterval = 10;
    scan.screen_dist = 5000; //mm
    vector<double> scandists;

    //arguments2
    string path(".");
    bool cvDepth(false), cvColor(false);
    bool pcd(false), pcd_rot(false), ply(false), ply_rot(false);
    for(int i=1;i<argc-2;i++){
        if(string(argv[i])=="-dist"){
            stringstream ss(argv[++i]);
            vector<double> inputs;
            double input;
            while(ss>>input) inputs.push_back(input);
            if(inputs.size()==1) scan.distance = inputs[0];
            else
                for(int n=0;n<inputs[2];n++)
                    scandists.push_back(inputs[0]+inputs[1]*n);
        }
        else if(string(argv[i])=="-rot")  scan.rotInterval = atoi(argv[++i]);
        else if(string(argv[i])=="-path") path = argv[++i];
        else if(string(argv[i])=="-cv"){
            if(string(argv[++i])=="depth") cvDepth = true;
            else if(string(argv[i])=="color") cvColor = true;
            else if(string(argv[i])=="both") {cvDepth = true; cvColor = true;}
        }
        else if(string(argv[i])=="-pcd"){
            pcd = true;
            if(atoi(argv[++i])==1) pcd_rot = true;
        }
        else if(string(argv[i])=="-ply"){
            ply = true;
            if(atoi(argv[++i])==1) ply_rot = true;
        }

        else PrintUsage();
    }
    if(scandists.size()==0) scandists.push_back(scan.distance);

    //read ply file
    vtkSmartPointer<vtkPolyData> data;
    vtkPLYReader* reader = vtkPLYReader::New ();
    reader->SetFileName (modelFileName.c_str ());
    reader->Update ();
    data = reader->GetOutput ();

    // Set the distance range
    double maxExt2(0);
    for(int i=0;i<data->GetNumberOfPoints();i++){
        double origin[3] = {0,0,0}, point[3];
        origin, data->GetPoint(i, point);
        double dist2 = vtkMath::Distance2BetweenPoints(origin, point);
        maxExt2 = dist2>maxExt2? dist2:maxExt2;
    }
    double maxExt  = sqrt(maxExt2);
    double minDist = scandists[0]-maxExt-100 < 0? 0:scandists[0]-maxExt-100;
    double maxDist = scandists.back() + maxExt + 100;
    double depthFactor = 2000./(maxDist-minDist) > 1? 1:2000./(maxDist-minDist);
    if(maxDist<2000) minDist = 0;
    cout<<"min. distance: "<<minDist<<endl;
    cout<<"max. distance: "<<maxDist<<endl;

    cout<<"Initializing..."<<flush;
    // Build a spatial locator for our dataset
    vtkSmartPointer<vtkCellLocator> tree = vtkSmartPointer<vtkCellLocator>::New ();
    tree->SetDataSet (data);
    tree->CacheCellBoundsOn ();
    tree->SetTolerance (0.0);
    tree->SetNumberOfCellsPerBucket (1);
    tree->AutomaticOn ();
    tree->BuildLocator ();
    tree->Update ();

    //read viewpoint ply file
    vtkSmartPointer<vtkPolyData> viewPointData;
    vtkPLYReader* viewPointReader = vtkPLYReader::New ();
    viewPointReader->SetFileName (viewFileName.c_str ());
    viewPointReader->Update ();
    viewPointData = viewPointReader->GetOutput ();

    //get viewpoints
    vector<double*> viewRays;
    for(vtkIdType i=0;i<viewPointData->GetNumberOfCells();i++){
        double a[3], b[3], c[3];
        viewPointData->GetPoint(viewPointData->GetCell(i)->GetPointId(0), a);
        viewPointData->GetPoint(viewPointData->GetCell(i)->GetPointId(1), b);
        viewPointData->GetPoint(viewPointData->GetCell(i)->GetPointId(2), c);

        double ab[3], ac[3];
        double* normal = new double[3];
        vtkMath::Subtract(a,b,ab);
        vtkMath::Subtract(a,c,ac);
        vtkMath::Cross(ac, ab, normal);
        vtkMath::Normalize(normal);
        viewRays.push_back(normal);
    }

    // Virtual camera parameters
    double viewray0[3] = {0.0, 0.0, 1.0};
    double up0[3]      = {0.0, 1.0, 0.0}; //normalize

    // Screen parameters
    double hor_len  = (scan.distance+scan.screen_dist) * 2. * tan(scan.fov_hor*0.5);
    double vert_len = (scan.distance+scan.screen_dist) * 2. * tan(scan.fov_vert*0.5);
    double hor_interval  = hor_len  / (double) scan.res_hor;
    double vert_interval = vert_len / (double) scan.res_vert;

    //create buffer for PLY/PCD print
    vector<double*> xyz;
    if(pcd||ply){
        xyz.reserve(scan.res_vert*scan.res_hor);
        for(int n=0;n<scan.res_vert*scan.res_hor;n++){
            xyz[n] = new double[3];
        }
    }

    //Initialize LINEMOD detector
    cv::Ptr<cv::linemod::Detector> detector;
    detector = cv::linemod::getDefaultLINEMOD();
    int num_modalities = (int)detector->getModalities().size();

    //OpenCV windows
    cv::namedWindow("depth");
    cv::namedWindow("color");
    cv::namedWindow("mask");
    Mouse::start("color");

    //main loop start
    if(path!=".") system(("mkdir "+path).c_str());
    ofstream log(path+"/labels.txt");
    log<<"maxDist "<<maxDist<<endl;
    log<<"minDist "<<minDist<<endl;
    log<<"model "<<modelFileName<<endl;
    log<<"viewpoint "<<viewFileName<<endl<<endl;
    vtkTransform* tr1 = vtkTransform::New ();
    vtkTransform* tr = vtkTransform::New ();
    vtkMatrix4x4* mat = vtkMatrix4x4::New();
    vtkMatrix4x4* mat1 = vtkMatrix4x4::New();

    int totScans = viewRays.size()*scandists.size()*360./(double)scan.rotInterval;
    //int totScans = viewRays.size()*scandists.size()*11;
    int failCount(0), scanCount(0), skipCount(0);
    for(int i=0;i<viewRays.size();i++){ //viewpoint loop

        log<<"viewpoint "<<i<<" "<<viewRays[i][0]<<" "<<viewRays[i][1]<<" "<<viewRays[i][2]<<endl<<endl;
        //set viewray and upVector
        double angle = vtkMath::AngleBetweenVectors(viewRays[i], viewray0); //radian
        double axis[3], up1[3];
        vtkMath::Cross(viewray0,viewRays[i], axis);
        tr1->Identity ();
        tr1->RotateWXYZ(angle*180/pi, axis); //requires angle in degree
        tr1->TransformPoint(up0, up1);
        tr1->Inverse();
        tr1->GetMatrix(mat1);
        vtkMath::Normalize(up1);
        bool skip(false);

        for(int deg = 0; deg<360; deg+=scan.rotInterval){ //rotation loop
            bool degSkip = false;

            //set up/right vectors
            double up[3], right[3];
            tr->Identity ();
            tr->RotateWXYZ((double)deg, viewRays[i]);
            tr->TransformPoint(up1, up);
            tr->Inverse();
            tr->GetMatrix(mat);
            vtkMatrix4x4::Multiply4x4(mat1,mat, mat);
            tr->Identity();
            tr->SetMatrix(mat);

            vtkMath::Normalize(up);
            vtkMath::Cross(viewRays[i], up, right);
            vtkMath::Normalize(right);

            for(int dist=0;dist<scandists.size();dist++){ //distance loop
                cout<<"\rScanning... "<<++scanCount<<"/"<<totScans-skipCount<<"          "<<flush;
                scan.distance = scandists[dist];
                mat->Element[2][3] = scan.distance;
                //Start!!
                double eye[3];
                eye[0] = -viewRays[i][0]*scan.distance;
                eye[1] = -viewRays[i][1]*scan.distance;
                eye[2] = -viewRays[i][2]*scan.distance;

                // Generate screen
                double screen_center[3];
                for(int n=0;n<3;n++) screen_center[n] = viewRays[i][n] * scan.screen_dist;
                double screen_leftup[3];
                for(int n=0;n<3;n++) screen_leftup[n] = screen_center[n]+up[n]*vert_len*0.5 - right[n]*hor_len*0.5;

                pcl::PointCloud<pcl::PointXYZ> screenCloud;
                for(int vert = 0 ; vert<scan.res_vert;vert++){
                    double left[3] = {screen_leftup[0] - up[0]*vert_interval*(vert+0.5),
                                      screen_leftup[1] - up[1]*vert_interval*(vert+0.5),
                                      screen_leftup[2] - up[2]*vert_interval*(vert+0.5)};
                    for(int hor = 0 ; hor<scan.res_hor;hor++){
                        double point[3] = {left[0]+right[0]*hor_interval*(hor+0.5),
                                           left[1]+right[1]*hor_interval*(hor+0.5),
                                           left[2]+right[2]*hor_interval*(hor+0.5)};
                        screenCloud.push_back(pcl::PointXYZ(point[0], point[1], point[2]));
                    }
                }

                // Scanning
                vector<int> rr, gg, bb;
                double p_coords[3], x[3], t, rgb[3];
                //create ::Mat for depth & color
                cv::Mat depth(scan.res_vert, scan.res_hor, CV_16U, cv::Scalar::all(0));
                cv::Mat color(scan.res_vert, scan.res_hor, CV_8UC3, cv::Scalar(0,0,0));
                uchar* color_data = color.data;
                cv::Mat mask(scan.res_vert, scan.res_hor, CV_8U, cv::Scalar::all(0));
                uchar* mask_data = mask.data;
                int subId;
                for(int vert = 0 ; vert<scan.res_vert;vert++){
                    for(int hor = 0 ; hor<scan.res_hor;hor++){
                        vtkIdType cellId;
                        double point[3];
                        point[0] = screenCloud.at(vert*scan.res_hor+hor).x;
                        point[1] = screenCloud.at(vert*scan.res_hor+hor).y;
                        point[2] = screenCloud.at(vert*scan.res_hor+hor).z;
                        if (tree->IntersectWithLine (eye, point, 0, t, x, p_coords, subId, cellId))
                        {
                            double ep[3];
                            vtkMath::Subtract(x, eye, ep);
                            depth.at<ushort>(vert,hor) = vtkMath::Dot(ep, viewRays[i]);
                            mask_data[vert*scan.res_hor+hor] = 255;
                            EstimateColor(data, cellId, x, rgb);
                            color_data[vert*scan.res_hor*3+hor*3+0] = floor(rgb[2]+0.5);
                            color_data[vert*scan.res_hor*3+hor*3+1] = floor(rgb[1]+0.5);
                            color_data[vert*scan.res_hor*3+hor*3+2] = floor(rgb[0]+0.5);
                             if(ply||pcd){
                                xyz[vert*scan.res_hor+hor][0]=x[0];
                                xyz[vert*scan.res_hor+hor][1]=x[1];
                                xyz[vert*scan.res_hor+hor][2]=x[2];
                                if(ply){
                                    rr.push_back(floor(rgb[0]+0.5));gg.push_back(floor(rgb[1]+0.5));bb.push_back(floor(rgb[2]+0.5));
                                }
                            }
                        }
                        else if(ply||pcd){
                            xyz[vert*scan.res_hor+hor][0]=point[0];
                            xyz[vert*scan.res_hor+hor][1]=point[1];
                            xyz[vert*scan.res_hor+hor][2]=point[2];
                            if(ply){
                                rr.push_back(0);gg.push_back(0);bb.push_back(0);
                            }
                        }
                    } // Horizontal
                } // Vertical

                //addTemplate (LINEMOD)
                depth = (depth-minDist)*depthFactor;
                std::vector<cv::Mat> sources;
                sources.push_back(color);
                sources.push_back(depth);
                cv::Rect bBox;
                int template_id = detector->addTemplate(sources, to_string(i), mask, &bBox);
                cv::Mat display = color;
                if(template_id>=0){
                    const std::vector<cv::linemod::Template>& templates = detector->getTemplates(to_string(i), template_id);
                    drawResponse(templates, num_modalities, display, cv::Point(bBox.x, bBox.y), detector->getT(0));
                    //log<<to_string(i)+" "+to_string(deg)<<" "<<"upVector: "<<up[0]<<" "<<up[1]<<" "<<up[2]<<endl;
                    log<<"template "<<template_id<<endl;
                    log<<"distance "<<scan.distance<<endl;
                    log<<"rotation "<<endl;
                    log<<*mat<<endl<<endl;
                }else failCount++;
                std::vector<cv::linemod::Match> matches;
                std::vector<cv::String> class_ids;
                std::vector<cv::Mat> quantized_images;
                detector->match(sources, 80., matches, class_ids, quantized_images);
                //visualizing
                cv::imshow("depth", quantized_images[1]);
                cv::putText(display, "view: "+to_string(i)+"/"+to_string(viewRays.size())+" (deg: "+to_string(deg)+" dist: "+to_string((int)scan.distance)+")"
                            , cv::Point(3,15), cv::FONT_HERSHEY_SIMPLEX, 0.5, cv::Scalar(255,255,255),0.5);
                cv::imshow("color", display);
                cv::imshow("mask", mask);

                //printing
                //** 1)OpenCV
                if(cvDepth){
                    std::string fName(path + "/" + to_string(i)+"_r"+to_string(deg)+"_d"+to_string((int)scan.distance)+".yml");
                    cv::FileStorage fs(fName, cv::FileStorage::WRITE);
                    fs << "depth" << depth; fs.release();
                }
                //** 2)PLY
                if(ply){
                    ofstream ofs(path+"/"+to_string(i)+"_r"+to_string(deg)+"_d"+to_string((int)scan.distance)+".ply");
                    ofs<<"ply"<<endl;
                    ofs<<"format ascii 1.0"<<endl;
                    ofs<<"comment exported in vitual_scanner"<<endl;
                    ofs<<"element vertex "<<scan.res_vert*scan.res_hor<<endl;
                    ofs<<"property float x"<<endl;
                    ofs<<"property float y"<<endl;
                    ofs<<"property float z"<<endl;
                    ofs<<"property uchar red"<<endl;
                    ofs<<"property uchar green"<<endl;
                    ofs<<"property uchar blue"<<endl;
                    ofs<<"end_header"<<endl;
                    for(int idx=0;idx<scan.res_vert*scan.res_hor;idx++){
                        if(ply_rot){
                            double xyz_rot[3];
                            tr->TransformPoint(xyz[idx],xyz_rot);
                            ofs<<xyz_rot[0]<<" "<<xyz_rot[1]<<" "<<xyz_rot[2]<<" "<<rr[idx]<<" "<<gg[idx]<<" "<<bb[idx]<<endl;
                        }
                        else
                            ofs<<xyz[idx][0]<<" "<<xyz[idx][1]<<" "<<xyz[idx][2]<<" "<<rr[idx]<<" "<<gg[idx]<<" "<<bb[idx]<<endl;
                    }
                    ofs.close();
                }

                //** 3)PCD
                if(pcd){
                    pcl::PointCloud<pcl::PointXYZ> cloud;
                    for(int idx=0;idx<scan.res_vert*scan.res_hor;idx++){
                        if(ply_rot){
                            double xyz_rot[3];
                            tr->TransformPoint(xyz[idx],xyz_rot);
                            cloud.push_back(PointXYZ(xyz_rot[0],xyz_rot[1],xyz_rot[2]));
                        }
                        else
                            cloud.push_back(PointXYZ(xyz[idx][0],xyz[idx][1],xyz[idx][2]));
                    }
                    pcl::PCDWriter writer;
                    writer.writeBinaryCompressed(path+"/"+to_string(i)+"_r"+to_string(deg)+"_d"+to_string((int)scan.distance)+".pcd", cloud);
                }
                char key = (char)cv::waitKey(10);
                if (key == 'q')
                    break;

                switch (key)
                {
                case 'h':
                    PrintKeyUsage();
                    break;
                case 'w':
                    writeLinemod(detector, path+"/"+"templates.yml");
                    cout<<endl<<"w -> printed the template in "+path+"/"+"templates.yml"<<endl;
                    break;
                case 'v':
                    skip = true;
                    skipCount += scandists.size()*(floor((360-deg)/scan.rotInterval)-1)+scandists.size()-dist-1;
                    cout<<endl<<"v -> to next view"<<endl;
                    break;
                case 'd':
                    degSkip = true;
                    skipCount += scandists.size()-dist-1;
                    cout<<endl<<"v -> to next view"<<endl;
                    break;
                case 'f':
                    cout<<endl<<"f -> # of failed templates: "<<failCount<<endl;
                    break;
                case 'p':
                    cout<<endl<<"p -> pause press any key "<<failCount<<endl;
                    cv::waitKey();
                    break;

                default:
                    ;
                }
                if(skip) break;
                if(degSkip) break;
            }if(skip) break;
        }
    }log.close();
    writeLinemod(detector, path+"/"+"templates.yml");
    cout<<endl<<"failed for "<<failCount<<" templates"<<endl;
    return 0;
}

void EstimateColor(vtkPolyData* poly, vtkIdType cellId, double* point, double* rgb){
    double a[3], b[3], c[3];
    poly->GetPoint(poly->GetCell(cellId)->GetPointId(0), a);
    poly->GetPoint(poly->GetCell(cellId)->GetPointId(1), b);
    poly->GetPoint(poly->GetCell(cellId)->GetPointId(2), c);
    double ab[3], ac[3], ba[3], bc[3];
    vtkMath::Subtract(b, a, ab);
    vtkMath::Subtract(c, a, ac);
    vtkMath::Subtract(a, b, ba);
    vtkMath::Subtract(c, b, bc);

    double n0[3], n[3], end[3], l[3], d, l0p0[3];
    //a (p0=b, l0=a, l=ap, n0=(baXbc), n=n0Xbc)
    vtkMath::Subtract(point, a, l); //l=ap
    vtkMath::Cross(ba, bc, n0);     //n0
    vtkMath::Cross(n0, bc, n);      //n
    vtkMath::Subtract(b, a, l0p0);  //l0p0=p0-l0
    d = vtkMath::Dot(l0p0, n)/vtkMath::Dot(l, n);
    vtkMath::MultiplyScalar(l, d);
    vtkMath::Add(a, l, end);
    double a_w = sqrt(vtkMath::Distance2BetweenPoints(point, end)/vtkMath::Distance2BetweenPoints(a, end));
    //b (p0=a, l0=b, l=bp, n0=(abXac), n=n0Xac)
    vtkMath::Subtract(point, b, l); //l=bp
    vtkMath::Cross(ab, ac, n0);     //n0
    vtkMath::Cross(n0, ac, n);      //n
    vtkMath::Subtract(a, b, l0p0);  //l0p0=p0-l0
    d = vtkMath::Dot(l0p0, n)/vtkMath::Dot(l, n);
    vtkMath::MultiplyScalar(l, d);
    vtkMath::Add(b, l, end);
    double b_w = sqrt(vtkMath::Distance2BetweenPoints(point, end)/vtkMath::Distance2BetweenPoints(b, end));
    //c (p0=a, l0=c, l=cp, n0=(baXbc), n=n0Xba)
    vtkMath::Subtract(point, c, l); //l=cp
    vtkMath::Cross(ba, bc, n0);     //n0
    vtkMath::Cross(n0, ba, n);      //n
    vtkMath::Subtract(a, c, l0p0);  //l0p0=p0-l0
    d = vtkMath::Dot(l0p0, n)/vtkMath::Dot(l, n);
    vtkMath::MultiplyScalar(l, d);
    vtkMath::Add(c, l, end);
    double c_w = sqrt(vtkMath::Distance2BetweenPoints(point, end)/vtkMath::Distance2BetweenPoints(c, end));

    double rgb_a[3] = {poly->GetAttributes(0)->GetArray(0)->GetComponent(poly->GetCell(cellId)->GetPointId(0),0),
                       poly->GetAttributes(0)->GetArray(0)->GetComponent(poly->GetCell(cellId)->GetPointId(0),1),
                       poly->GetAttributes(0)->GetArray(0)->GetComponent(poly->GetCell(cellId)->GetPointId(0),2)};
    double rgb_b[3] = {poly->GetAttributes(0)->GetArray(0)->GetComponent(poly->GetCell(cellId)->GetPointId(1),0),
                       poly->GetAttributes(0)->GetArray(0)->GetComponent(poly->GetCell(cellId)->GetPointId(1),1),
                       poly->GetAttributes(0)->GetArray(0)->GetComponent(poly->GetCell(cellId)->GetPointId(1),2)};
    double rgb_c[3] = {poly->GetAttributes(0)->GetArray(0)->GetComponent(poly->GetCell(cellId)->GetPointId(2),0),
                       poly->GetAttributes(0)->GetArray(0)->GetComponent(poly->GetCell(cellId)->GetPointId(2),1),
                       poly->GetAttributes(0)->GetArray(0)->GetComponent(poly->GetCell(cellId)->GetPointId(2),2)};

    vtkMath::MultiplyScalar(rgb_a, a_w);
    vtkMath::MultiplyScalar(rgb_b, b_w);
    vtkMath::MultiplyScalar(rgb_c, c_w);
    vtkMath::Add(rgb_a, rgb_b, rgb);
    vtkMath::Add(rgb, rgb_c, rgb);
}

void drawResponse(const std::vector<cv::linemod::Template>& templates,
    int num_modalities, cv::Mat& dst, cv::Point offset, int T)
{
    static const cv::Scalar COLORS[5] = { CV_RGB(0, 0, 255),
                                          CV_RGB(0, 255, 0),
                                          CV_RGB(255, 255, 0),
                                          CV_RGB(255, 140, 0),
                                          CV_RGB(255, 0, 0) };

    for (int m = 0; m < num_modalities; ++m)
    {
        // NOTE: Original demo recalculated max response for each feature in the TxT
        // box around it and chose the display color based on that response. Here
        // the display color just depends on the modality.
        cv::Scalar color = COLORS[m];

        for (int i = 0; i < (int)templates[m].features.size(); ++i)
        {
            cv::linemod::Feature f = templates[m].features[i];
            cv::Point pt(f.x + offset.x, f.y + offset.y);
            cv::circle(dst, pt, T / 2, color);
        }
    }
}

static void writeLinemod(const cv::Ptr<cv::linemod::Detector>& detector, const std::string& filename)
{
    cv::FileStorage fs(filename, cv::FileStorage::WRITE);
    detector->write(fs);

    std::vector<cv::String> ids = detector->classIds();
    fs << "classes" << "[";
    for (int i = 0; i < (int)ids.size(); ++i)
    {
        fs << "{";
        detector->writeClass(ids[i], fs);
        fs << "}"; // current class
    }
    fs << "]"; // classes
}
