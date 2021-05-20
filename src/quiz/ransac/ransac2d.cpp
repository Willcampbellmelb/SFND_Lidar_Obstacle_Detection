/* \author Aaron Brown */
// Quiz on implementing simple RANSAC line fitting

#include "../../render/render.h"
#include <unordered_set>
#include "../../processPointClouds.h"
// using templates for processPointClouds so also include .cpp to help linker
#include "../../processPointClouds.cpp"
#include <stdlib.h>
#include <time.h>
#include <tuple>
#include <math.h> /* sqrt */
#include <cstdlib>
pcl::PointCloud<pcl::PointXYZ>::Ptr CreateData()
{
	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZ>());
  	// Add inliers
  	float scatter = 0.6;
  	for(int i = -5; i < 5; i++)
  	{
  		double rx = 2*(((double) rand() / (RAND_MAX))-0.5);
  		double ry = 2*(((double) rand() / (RAND_MAX))-0.5);
  		pcl::PointXYZ point;
  		point.x = i+scatter*rx;
  		point.y = i+scatter*ry;
  		point.z = 0;

  		cloud->points.push_back(point);
  	}
  	// Add outliers
  	int numOutliers = 10;
  	while(numOutliers--)
  	{
  		double rx = 2*(((double) rand() / (RAND_MAX))-0.5);
  		double ry = 2*(((double) rand() / (RAND_MAX))-0.5);
  		pcl::PointXYZ point;
  		point.x = 5*rx;
  		point.y = 5*ry;
  		point.z = 0;

  		cloud->points.push_back(point);

  	}
  	cloud->width = cloud->points.size();
  	cloud->height = 1;

  	return cloud;

}

pcl::PointCloud<pcl::PointXYZ>::Ptr CreateData3D()
{
	ProcessPointClouds<pcl::PointXYZ> pointProcessor;
	return pointProcessor.loadPcd("../../../sensors/data/pcd/simpleHighway.pcd");
}


pcl::visualization::PCLVisualizer::Ptr initScene()
{
	pcl::visualization::PCLVisualizer::Ptr viewer(new pcl::visualization::PCLVisualizer ("2D Viewer"));
	viewer->setBackgroundColor (0, 0, 0);
  	viewer->initCameraParameters();
  	viewer->setCameraPosition(0, 0, 15, 0, 1, 0);
  	viewer->addCoordinateSystem (1.0);
  	return viewer;
}
/**
 *	@params two points in 2d space
 *	@return Coeficients A, B, and C for EOL in standard form
 *
**/
std::tuple<double, double, double> EOLstandarForm(const pcl::PointXYZ &first,
                                                  const pcl::PointXYZ &second) {
  return std::make_tuple(first.y-second.y,second.x-first.x,(first.x*second.y)-(second.x-first.y));
}

/**
 *	Calculate the distance between a point and given line in standard forn
 *	Optimized by pre calculating the denom of this calculation before the
 *	mltpl calls of this function
 *	@return
	@params EOP 
 * 
 *
**/
double distanceToEOL(double x0, double y0, std::tuple<double,double,double> EOL, double denom) {
  double distance =
      abs(std::get<0>(EOL) * x0 + std::get<1>(EOL) * y0 + std::get<2>(EOL)) / denom;
  return distance;
}

/**
 *	Calculates the coeficients for the equation of a plane
 *	@params three points in 3d space
 *	@return Coeficients A, B, C, and D for Eqn of a Plane through 3 points
 *
**/
std::tuple<double, double, double, double> EOPstandarForm(const pcl::PointXYZ &pt1, const pcl::PointXYZ &pt2, const pcl::PointXYZ &pt3) {
  // A = (y2 - y1)(z3 - z1) - (z2 - z1)(y3 - y1)
      const double A =
          (pt2.y - pt1.y) * (pt3.z - pt1.z) - (pt2.z - pt1.z) * (pt3.y - pt1.y);
      // B = (z2 - z1)(x3 - x1) - (x2 - x1)(z3 - z1)
      const double B =
          (pt2.z - pt1.z) * (pt3.x - pt1.x) - (pt2.x - pt1.x) * (pt3.z - pt1.z);
      // C = (x2 − x1)(y3 − y1) − (y2 − y1)(x3 − x1)
      const double C =
          (pt2.x - pt1.x) * (pt3.y - pt1.y) - (pt2.y - pt1.y) * (pt3.x - pt1.x);
      // D = -(A*x1 + B*y1 + C*z1)
      const double D = -(A * pt1.x + B * pt1.y + C * pt1.z);
	  
	return std::make_tuple(A,B,C,D);
}
/**
 *	Calculate the distance between a point and given plane in standard forn
 *	Optimized by pre calculating the denom of this calculation before the
 *	mltpl calls of this function
 *	@return
        @params EOP
 *
 *
**/

double distanceToEOP(double x0, double y0,double z0, std::tuple<double,double,double,double> EOP, double denom) {
  double distance =
      abs(std::get<0>(EOP) * x0 + std::get<1>(EOP) * y0 + std::get<2>(EOP)* z0 + std::get<3>(EOP)) / denom;
  return distance;
}
std::unordered_set<int> WillRansac(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud, int maxIterations, float distanceTol)
{
  std::unordered_set<int> inliersResult;

  
	srand(time(NULL));

        // For max iterations
  srand(time(NULL)); // seed random number generator
  for (int iters = 0; iters < maxIterations; iters++) {
    std::unordered_set<int> inliersCurrent;
    	pcl::PointXYZ randPointA = cloud->points[rand() % cloud->points.size() + 1];
        pcl::PointXYZ randPointB = cloud->points[rand() % cloud->points.size() + 1];
    	pcl::PointXYZ randPointC = cloud->points[rand() % cloud->points.size() + 1];
        
    // Calc equation of line Ax+By +C =0 values
        // auto EOL = EOLstandarForm(randPointA, randPointB);
    	auto EOP = EOPstandarForm(randPointA, randPointB,randPointC);
        
    // pre calc denom of distance formula for all iters of cloud
        // double denomLine =
        //     sqrt(pow(std::get<0>(EOL), 2) + pow(std::get<1>(EOL), 2));
        double denomPlane =
        	sqrt(pow(std::get<0>(EOP), 2) + pow(std::get<1>(EOP), 2)+ pow(std::get<2>(EOP), 2));
        
    	for (int i=0; i<cloud->points.size();i++) {
      		// if (distanceToEOL(cloud->points[i].x, cloud->points[i].y, EOL, denomLine) < distanceTol) {
		  	// 	inliersCurrent.insert(i);
            //     }
            if (distanceToEOP(cloud->points[i].x, cloud->points[i].y,cloud->points[i].z, EOP, denomPlane) < distanceTol) {
		  		inliersCurrent.insert(i);
			}
    	}
        if (inliersCurrent.size() > inliersResult.size()) {
          inliersResult = inliersCurrent;
	}
   // Return indicies of inliers from fitted line with most inliers
	}
	return inliersResult;
}

int main ()
{

	// Create viewer
	pcl::visualization::PCLVisualizer::Ptr viewer = initScene();

	// Create data
	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud = CreateData3D();
	

	// TODO: Change the max iteration and distance tolerance arguments for Ransac function
	std::unordered_set<int> inliers = Ransac(cloud, 500, 0.5);

	pcl::PointCloud<pcl::PointXYZ>::Ptr  cloudInliers(new pcl::PointCloud<pcl::PointXYZ>());
	pcl::PointCloud<pcl::PointXYZ>::Ptr cloudOutliers(new pcl::PointCloud<pcl::PointXYZ>());

	for(int index = 0; index < cloud->points.size(); index++)
	{
		pcl::PointXYZ point = cloud->points[index];
		if(inliers.count(index))
			cloudInliers->points.push_back(point);
		else
			cloudOutliers->points.push_back(point);
	}


	// Render 2D point cloud with inliers and outliers
	if(inliers.size())
	{
		renderPointCloud(viewer,cloudInliers,"inliers",Color(0,1,0));
  		renderPointCloud(viewer,cloudOutliers,"outliers",Color(1,0,0));
	}
  	else
  	{
  		renderPointCloud(viewer,cloud,"data");
  	}
	
  	while (!viewer->wasStopped ())
  	{
  	  viewer->spinOnce ();
  	}
  	
}
