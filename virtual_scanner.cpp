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

#include "functions.cpp"

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
    scan.distance = 2400; //mm
    scan.rotInterval = 10;
    scan.screen_dist = 5000; //mm
    vector<double> scandists;

    //arguments2
    string path(".");
    bool cvDepth(false), cvColor(false);
    bool pcd(false), pcd_rot(false), ply(false), ply_rot(false);
    for(int i=1;i<argc-2;i++){
		if (string(argv[i]) == "-dist") {
			vector<double> inputs;
			for (int j=0; j<3;j ++) inputs.push_back(atoi(argv[++i]));
			double input;
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
    vtkPLYReader* reader = vtkPLYReader::New();
    reader->SetFileName(modelFileName.c_str());
    reader->Update();
    data = reader->GetOutput(); // data: 3D scan model

    // Set the depth factor by calculating model max/min distance
    double maxExt2(0);
    for(int i=0;i<data->GetNumberOfPoints();i++){
        double origin[3] = {0,0,0}, point[3];
        data->GetPoint(i, point);
        double dist2 = vtkMath::Distance2BetweenPoints(origin, point); // distance^2 between orgin and vertex
        maxExt2 = dist2>maxExt2? dist2:maxExt2;
    }
    double maxExt  = sqrt(maxExt2);
    double minDist = scandists[0]-maxExt-100 < 0? 0:scandists[0]-maxExt-100;
    double maxDist = scandists.back() + maxExt + 100;
    double depthFactor = 2000./(maxDist-minDist) > 1? 1:2000./(maxDist-minDist);

    if(maxDist<2000) minDist = 0;
    cout<<"min. distance: "<<minDist<<endl;
    cout<<"max. distance: "<<maxDist<<endl;
    // ------------------------------------

    cout<<"Initializing..."<<flush;
    // Build a spatial locator for our dataset
    vtkSmartPointer<vtkCellLocator> tree = vtkSmartPointer<vtkCellLocator>::New();
    tree->SetDataSet(data);
    tree->CacheCellBoundsOn();
    tree->SetTolerance(0.0);
    tree->SetNumberOfCellsPerBucket(1);
    tree->AutomaticOn();
    tree->BuildLocator();
    tree->Update();

    // Read viewpoint ply file
    vtkSmartPointer<vtkPolyData> viewPointData;
    vtkPLYReader* viewPointReader = vtkPLYReader::New();
    viewPointReader->SetFileName (viewFileName.c_str());
    viewPointReader->Update();
    viewPointData = viewPointReader->GetOutput(); // viewpoint model data

    // Get viewpoints using polyhedral sphere facets

    vector<double*> viewRays; // camera viewpoints (inside normals)
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
    double up0[3]      = {0.0, -1.0, 0.0}; //normalize

    // Screen parameters
    double hor_len  = (scan.distance+scan.screen_dist) * 2. * tan(scan.fov_hor*0.5);
    double vert_len = (scan.distance+scan.screen_dist) * 2. * tan(scan.fov_vert*0.5);
    double hor_interval  = hor_len  / (double) scan.res_hor;
    double vert_interval = vert_len / (double) scan.res_vert;

    // Create buffer for PLY/PCD print
    vector<double*> xyz;
    if(pcd||ply){
        xyz.reserve(scan.res_vert*scan.res_hor);
        for(int n=0;n<scan.res_vert*scan.res_hor;n++){
            xyz[n] = new double[3];
        }
    }

    // Initialize LINEMOD detector
    cv::Ptr<cv::linemod::Detector> detector;
    detector = cv::linemod::getDefaultLINEMOD();
    int num_modalities = (int)detector->getModalities().size();

    // OpenCV windows
    cv::namedWindow("depth");
    cv::namedWindow("color");
    cv::namedWindow("mask");
    Mouse::start("color");

    

    //=====================================================================================
    // Main loop start ====================================================================
    //=====================================================================================
    PrintKeyUsage();

    if(path!=".") system(("mkdir "+path).c_str());
    ofstream log(path+"/labels.txt");
    log<<"maxDist "<<maxDist<<endl;
    log<<"minDist "<<minDist<<endl;
    log<<"model "<<modelFileName<<endl;
    log<<"viewpoint "<<viewFileName<<endl<<endl;
    vtkTransform* tr1 = vtkTransform::New();  // Transform for the viewpoints
    vtkTransform* tr = vtkTransform::New();   // Transform for the angles
    vtkMatrix4x4* mat = vtkMatrix4x4::New();  // Affine Transform Matrix for the angles
    vtkMatrix4x4* mat1 = vtkMatrix4x4::New(); // Affine Transform Matrix for the viewpoints

    int totScans = viewRays.size()*scandists.size()*360./(double)scan.rotInterval;
    //int totScans = viewRays.size()*scandists.size()*11;
    int failCount(0), scanCount(0), skipCount(0);

    for(int i=0;i<viewRays.size();i++){ // 1. Viewpoint loop ==============================
    // for(int i=118;i<180;i++){ // 1. Viewpoint loop ==============================

        log<<"viewpoint "<<i<<" "<<viewRays[i][0]<<" "<<viewRays[i][1]<<" "<<viewRays[i][2]<<endl<<endl;
        // Set viewRay and upVector
        double angle = vtkMath::AngleBetweenVectors(viewRays[i], viewray0); // radian
        double axis[3], up1[3];
        vtkMath::Cross(viewray0,viewRays[i], axis);
        tr1->Identity ();
        tr1->RotateWXYZ(angle*180/pi, axis); // Requires angle in degree
        tr1->TransformPoint(up0, up1);       // Transform upvector wathcing the object
        tr1->Inverse();                      // For fixed object viewpoint 
        tr1->GetMatrix(mat1);                // inverse R1
        vtkMath::Normalize(up1);             // Normalize upVector
        bool skip(false);

        for(int deg = 0; deg<360; deg+=scan.rotInterval){ // 2. Rotation loop =============
        // for(int deg = 0; deg<60; deg+=scan.rotInterval){ // 2. Rotation loop =============
            bool degSkip = false;

            //set up/right vectors
            double up[3], right[3];
            tr->Identity();
            tr->RotateWXYZ((double)deg, viewRays[i]);
            tr->TransformPoint(up1, up); // Transform 'upvector for the axis of viewpoint
            tr->Inverse();               // For fixed object rotation
            tr->GetMatrix(mat);          // inverse R2
            vtkMatrix4x4::Multiply4x4(mat1,mat, mat);  // R = inv R1 * inv R2
            // tr->Identity();
            tr->SetMatrix(mat);

            vtkMath::Normalize(up);      // Normalize upVector
            vtkMath::Cross(viewRays[i], up, right);
            // vtkMath::Cross(up, viewRays[i], right);
            vtkMath::Normalize(right);

            for(int dist=0;dist<scandists.size();dist++){ // 3. Distance loop =============
                cout<<"\rScanning... "<<++scanCount<<"/"<<totScans-skipCount<<"          "<<flush;
                scan.distance = scandists[dist];
                mat->Element[2][3] = scan.distance;
                // Start!!
                double eye[3]; // object to facets (outside normal)
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
                // create ::Mat for depth & color (1280 x 720)
                cv::Mat depth(scan.res_vert, scan.res_hor, CV_16U, cv::Scalar::all(0));
                cv::Mat color(scan.res_vert, scan.res_hor, CV_8UC3, cv::Scalar(0,0,0));
                uchar* color_data = color.data;
                cv::Mat mask(scan.res_vert, scan.res_hor, CV_8U, cv::Scalar::all(0));
                uchar* mask_data = mask.data;

                int leftupVert(scan.res_vert);
                int leftupHor(scan.res_hor);
                vector<tuple<int,int,int>> rescheck;

                int subId;
                // Generate point intersected with line
                for(int vert = 0 ; vert<scan.res_vert;vert++){
                    for(int hor = 0 ; hor<scan.res_hor;hor++){
                        vtkIdType cellId;
                        double point[3];
                        // screen data
                        point[0] = screenCloud.at(vert*scan.res_hor+hor).x;
                        point[1] = screenCloud.at(vert*scan.res_hor+hor).y;
                        point[2] = screenCloud.at(vert*scan.res_hor+hor).z;
                        if (tree->IntersectWithLine (eye, point, 0, t, x, p_coords, subId, cellId))
                        {
                            double ep[3]; // eye to point
                            vtkMath::Subtract(x, eye, ep); // x: intersect point
                            depth.at<ushort>(vert,hor) = vtkMath::Dot(ep, viewRays[i]); // depth calculation with projection
                            mask_data[vert*scan.res_hor+hor] = 255;

                            // RGB color calculated with barycentric coordinate
                            EstimateColor(data, cellId, x, rgb);
                            color_data[vert*scan.res_hor*3+hor*3+0] = floor(rgb[2]+0.5);
                            color_data[vert*scan.res_hor*3+hor*3+1] = floor(rgb[1]+0.5);
                            color_data[vert*scan.res_hor*3+hor*3+2] = floor(rgb[0]+0.5);
                            // --------------------------------------------------------

                            // set leftup edgepoint
                            if (vert < leftupVert) leftupVert = vert;
                            if (hor  < leftupHor ) leftupHor  = hor;

                            // rescheck.push_back(make_tuple(hor,vert,10));

                             if(ply||pcd){ // set intersect point
                                xyz[vert*scan.res_hor+hor][0]=x[0];
                                xyz[vert*scan.res_hor+hor][1]=x[1];
                                xyz[vert*scan.res_hor+hor][2]=x[2];
                                if(ply){
                                    rr.push_back(floor(rgb[0]+0.5));gg.push_back(floor(rgb[1]+0.5));bb.push_back(floor(rgb[2]+0.5));
                                }
                            }
                        }
                        else if(ply||pcd){ // set screen point
                            xyz[vert*scan.res_hor+hor][0]=point[0];
                            xyz[vert*scan.res_hor+hor][1]=point[1];
                            xyz[vert*scan.res_hor+hor][2]=point[2];
                            if(ply){
                                rr.push_back(0);gg.push_back(0);bb.push_back(0);
                            }
                        }
                        rescheck.push_back(make_tuple(hor,vert,0));
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
                    log<<"position "<< scan.res_hor*0.5-leftupHor << " " << scan.res_vert*0.5-leftupVert << endl;
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
                // **
                // std::string boxName(path + "/" + to_string(i)+"_r"+to_string(deg)+"_d"+to_string((int)scan.distance)+"_bbox.obj");
                // ofstream ofs_box(boxName);
                // ofs_box << "v " << leftupHor << " " << leftupVert << " 10" << endl;

                // for(int idx=0;idx<scan.res_vert*scan.res_hor;idx++){
                //     if (get<2>(rescheck[idx]) == 0) continue;
                //     ofs_box << "v " << get<0>(rescheck[idx]) << " "
                //                     << get<1>(rescheck[idx]) << " " 
                //                     << get<2>(rescheck[idx]) << endl;
                // }


                //** 1)OpenCV
                if(cvDepth){
                    std::string fName(path + "/" + to_string(i)+"_r"+to_string(deg)+"_d"+to_string((int)scan.distance)+".yml");
                    cv::FileStorage fs(fName, cv::FileStorage::WRITE);
                    fs << "depth" << depth; fs.release();
                }
                //** 2)PLY
                if(ply){
                    int leftup(leftupVert*scan.res_hor+leftupHor);
                    int ctr(scan.res_vert*0.5*scan.res_hor + scan.res_hor*0.5);
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
                    if (ply_rot) {
                        double xyz_rot[3], xyz_rot2[3];
                        tr->TransformPoint(xyz[leftup],xyz_rot);
                        tr->TransformPoint(xyz[ctr],xyz_rot2);
                        ofs<<xyz_rot[0]<<" "<<xyz_rot[1]<<" "<<xyz_rot[2]<<" "<<rr[leftup]<<" "<<gg[leftup]<<" "<<bb[leftup]<<endl;
                        ofs<<xyz_rot2[0]<<" "<<xyz_rot2[1]<<" "<<xyz_rot2[2]<<" "<<rr[leftup]<<" "<<gg[leftup]<<" "<<bb[leftup]<<endl;
                    } else {
                        ofs<<xyz[leftup][0]<<" "<<xyz[leftup][1]<<" "<<xyz[leftup][2]<<" "<<rr[leftup]<<" "<<gg[leftup]<<" "<<bb[leftup]<<endl;
                        ofs<<xyz[ctr][0]<<" "<<xyz[ctr][1]<<" "<<xyz[ctr][2]<<" "<<rr[leftup]<<" "<<gg[leftup]<<" "<<bb[leftup]<<endl;
                    }
                    
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

