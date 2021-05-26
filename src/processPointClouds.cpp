// PCL lib Functions for processing point clouds

#include "processPointClouds.h"
#define CUSTOMRANSAC


// constructor:
template <typename PointT> ProcessPointClouds<PointT>::ProcessPointClouds() {}

// de-constructor:
template <typename PointT> ProcessPointClouds<PointT>::~ProcessPointClouds() {}

template <typename PointT>
void ProcessPointClouds<PointT>::numPoints(typename pcl::PointCloud<PointT>::Ptr cloud) {
  std::cout << cloud->points.size() << std::endl;
}

template <typename PointT>
typename pcl::PointCloud<PointT>::Ptr ProcessPointClouds<PointT>::FilterCloud(typename pcl::PointCloud<PointT>::Ptr cloud, float filterRes, Eigen::Vector4f minPoint, Eigen::Vector4f maxPoint) {

  // Time segmentation process
  auto startTime = std::chrono::steady_clock::now();

  // TODO:: Fill in the function to do voxel grid point reduction and region
  // based filtering

  auto endTime = std::chrono::steady_clock::now();
  auto elapsedTime = std::chrono::duration_cast<std::chrono::milliseconds>(
      endTime - startTime);
  std::cout << "filtering took " << elapsedTime.count() << " milliseconds"
            << std::endl;

  return cloud;
}
/**
 *
 *
 *@return Point cloud ptr pair to obstacles and plane point cloud(road)
 **/
template <typename PointT>
std::pair<typename pcl::PointCloud<PointT>::Ptr, typename pcl::PointCloud<PointT>::Ptr> ProcessPointClouds<PointT>::SeparateClouds(pcl::PointIndices::Ptr inliers,typename pcl::PointCloud<PointT>::Ptr cloud) {

  // Create Point cloud ptr objects for extracting to
  typename pcl::PointCloud<PointT>::Ptr obstCloud(new pcl::PointCloud<PointT>), planeCloud(new pcl::PointCloud<PointT>);

  // the segmented plane road points ARE INLIERS
  for(auto index:inliers->indices){
    planeCloud->points.push_back(cloud->points[index]);

  }
  // Use PCL extract to create obstacles from the inliers
  pcl::ExtractIndices<PointT>extract;
  // Extract the inliers
  extract.setInputCloud(cloud);
  extract.setIndices(inliers);
  extract.setNegative(false);
  extract.filter(*planeCloud);
   // Extract the outliers
  extract.setNegative(true);  // True means points belonging to outliers remain
      // outliers are the points not belonging to the road
  extract.filter(*obstCloud);

  std::pair<typename pcl::PointCloud<PointT>::Ptr,typename pcl::PointCloud<PointT>::Ptr> segResult(obstCloud, planeCloud);
  return segResult;
}

/**
 *
 *
 *@return Point cloud ptr pair to obstacles and plane point cloud
 **/
template <typename PointT>
std::pair<typename pcl::PointCloud<PointT>::Ptr,typename pcl::PointCloud<PointT>::Ptr> ProcessPointClouds<PointT>::SegmentPlane(typename pcl::PointCloud<PointT>::Ptr cloud, int maxIterations,float distanceThreshold) {
  // Time segmentation process
  auto startTime = std::chrono::steady_clock::now();
  // Create the segmentation object
  pcl::SACSegmentation<PointT> seg;
  pcl::PointIndices::Ptr inliers(new pcl::PointIndices());
  pcl::ModelCoefficients::Ptr coefficients(new pcl::ModelCoefficients());

  // Optional
  seg.setOptimizeCoefficients(true);
  // Mandatory
  seg.setModelType(pcl::SACMODEL_PLANE);
  seg.setMethodType(pcl::SAC_RANSAC);
  seg.setMaxIterations(maxIterations);
  seg.setDistanceThreshold(distanceThreshold);

  // Segment the largest planar component from the remaining cloud
  seg.setInputCloud(cloud);
  seg.segment(*inliers, *coefficients);
  if (inliers->indices.size() == 0) {
  std::cerr << "Could not estimate a planar model for the given dataset."
              << std::endl;
  }


  auto segResult = SeparateClouds(inliers, cloud);
  auto endTime = std::chrono::steady_clock::now();
  auto elapsedTime = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);
  std::cout << "plane segmentation took " << elapsedTime.count() << " milliseconds" << std::endl;

  return segResult;
}
//** my segment plane function - should be stand in for the one above **/

template <typename PointT>
std::vector<typename pcl::PointCloud<PointT>::Ptr> ProcessPointClouds<PointT>::Clustering(typename pcl::PointCloud<PointT>::Ptr cloud, float clusterTolerance,int minSize, int maxSize) {

  // Time clustering process
  auto startTime = std::chrono::steady_clock::now();

  std::vector<typename pcl::PointCloud<PointT>::Ptr> clusters;

  typename pcl::search::KdTree<PointT>::Ptr tree (new pcl::search::KdTree<PointT>);
  tree->setInputCloud (cloud);

  std::vector<pcl::PointIndices> clusterIndices;
  pcl::EuclideanClusterExtraction<pcl::PointXYZ> ec;
  ec.setClusterTolerance (clusterTolerance); // 2cm
  ec.setMinClusterSize (minSize);
  ec.setMaxClusterSize (maxSize);
  ec.setSearchMethod (tree);
  ec.setInputCloud (cloud);
  ec.extract(clusterIndices);

  for (pcl::PointIndices currentIndices : clusterIndices) {
    typename pcl::PointCloud<PointT>::Ptr cloudCluster (new pcl::PointCloud<PointT>);
    for (int index : currentIndices.indices) {
      cloudCluster->points.push_back(cloud->points[index]);
    }
    cloudCluster->width = cloudCluster->size ();
    cloudCluster->height = 1;
    cloudCluster->is_dense = true;
    clusters.push_back(cloudCluster);
  }
  //
  auto endTime = std::chrono::steady_clock::now();
  auto elapsedTime = std::chrono::duration_cast<std::chrono::milliseconds>(
      endTime - startTime);
  std::cout << "clustering took " << elapsedTime.count() << " milliseconds and found " << clusters.size() << " clusters" << std::endl;

  return clusters;
}

template <typename PointT>
Box ProcessPointClouds<PointT>::BoundingBox(
    typename pcl::PointCloud<PointT>::Ptr cluster) {

  // Find bounding box for one of the clusters
  PointT minPoint, maxPoint;
  pcl::getMinMax3D(*cluster, minPoint, maxPoint);

  Box box;
  box.x_min = minPoint.x;
  box.y_min = minPoint.y;
  box.z_min = minPoint.z;
  box.x_max = maxPoint.x;
  box.y_max = maxPoint.y;
  box.z_max = maxPoint.z;

  return box;
}

template <typename PointT>
void ProcessPointClouds<PointT>::savePcd(
    typename pcl::PointCloud<PointT>::Ptr cloud, std::string file) {
  pcl::io::savePCDFileASCII(file, *cloud);
  std::cerr << "Saved " << cloud->points.size() << " data points to " + file
            << std::endl;
}

template <typename PointT>
typename pcl::PointCloud<PointT>::Ptr ProcessPointClouds<PointT>::loadPcd(std::string file) {

  typename pcl::PointCloud<PointT>::Ptr cloud(new pcl::PointCloud<PointT>);

  if (pcl::io::loadPCDFile<PointT>(file, *cloud) == -1) //* load the file
  {
    PCL_ERROR("Couldn't read file \n");
  }
  std::cerr << "Loaded " << cloud->points.size() << " data points from " + file
            << std::endl;

  return cloud;
}

template <typename PointT>
std::vector<boost::filesystem::path> ProcessPointClouds<PointT>::streamPcd(std::string dataPath) {
  
  std::vector<boost::filesystem::path> paths(boost::filesystem::directory_iterator{dataPath},boost::filesystem::directory_iterator{});

  // sort files in accending order so playback is chronological
  sort(paths.begin(), paths.end());

  return paths;
}
