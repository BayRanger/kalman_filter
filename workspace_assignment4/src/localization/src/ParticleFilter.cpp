#include "localization/ParticleFilter.h"
#include "localization/Util.h"

#include "tf/tf.h"

using namespace std;

ParticleFilter::ParticleFilter(int numberOfParticles) {
	this->numberOfParticles = numberOfParticles;

	// initialize particles
	for (int i = 0; i < numberOfParticles; i++) {
		this->particleSet.push_back(new Particle());
	}

	// this variable holds the estimated robot pose
	this->bestHypothesis = new Particle();

	// at each correction step of the filter only the laserSkip-th beam of a scan should be integrated
	this->laserSkip = 5;

	// distance map used for computing the likelihood field
	this->distMap = NULL;
}

ParticleFilter::~ParticleFilter() {
	// delete particles
	for (int i = 0; i < numberOfParticles; i++) {
		Particle* p = this->particleSet[i];
		delete p;
	}

	this->particleSet.clear();

	if (this->likelihoodField)
		delete[] this->likelihoodField;

	delete this->bestHypothesis;

	if (this->distMap)
		delete[] this->distMap;
}

int ParticleFilter::getNumberOfParticles() {
	return this->numberOfParticles;
}

std::vector<Particle*>* ParticleFilter::getParticleSet() {
	return &(this->particleSet);
}

void ParticleFilter::initParticlesUniform() {
    //get map properties
    int mapWidth, mapHeight;
    double mapResolution;
    this->getLikelihoodField(mapWidth, mapHeight,mapResolution);
	double weight =1.0/(numberOfParticles);
	for (int i = 0; i < numberOfParticles; i++) {
		double x_val = Util::uniformRandom(0,mapWidth * mapResolution);
		double y_val = Util::uniformRandom(0,mapHeight * mapResolution);
		double theta_val = Util::uniformRandom(-M_PI, M_PI);
		*(this->particleSet[i]) = Particle(x_val,y_val,theta_val,weight);
	}


}

void ParticleFilter::initParticlesGaussian(double mean_x, double mean_y,
		double mean_theta, double std_xx, double std_yy, double std_tt) {
	// TODO: here comes your code
	double weight =1.0/(numberOfParticles);
	for (int i = 0; i < numberOfParticles; i++) {
		//Particle* p;
	double x_val = Util::gaussianRandom(mean_x, std_xx);
	double y_val = Util::gaussianRandom(mean_y, std_yy);
	double theta_val = Util::gaussianRandom(mean_theta, std_tt);
	*(this->particleSet[i]) = Particle(x_val,y_val,theta_val,weight);
	}

}

/**
 *  Initializes the likelihood field as our sensor model.
 */
void ParticleFilter::setMeasurementModelLikelihoodField(
		const nav_msgs::OccupancyGrid& map, double zRand, double sigmaHit) {
	ROS_INFO("Creating likelihood field for laser range finder...");

	// create the likelihood field - with the same discretization as the occupancy grid map
	this->likelihoodField = new double[map.info.height * map.info.width];
	this->likelihoodFieldWidth = map.info.width;
	this->likelihoodFieldHeight = map.info.height;
	this->likelihoodFieldResolution = map.info.resolution;

    // calculates the distance map and stores it in member variable 'distMap'
	// for every map position it contains the distance to the nearest occupied cell.

    // Here you have to create your likelihood field
	// HINT0: sigmaHit is given in meters. You have to take into account the resolution of the likelihood field to apply it.
	// HINT1: You will need the distance map computed 3 lines above
	// HINT2: You can visualize it in the map_view when clicking on "show likelihood field" and "publish all".
	// HINT3: Storing probabilities in each cell between 0.0 and 1.0 might lead to round-off errors, therefore it is
	// good practice to convert the probabilities into log-space, i.e. storing log(p(x,y)) in each cell. As a further
	// advantage you can simply add the log-values in your sensor model, when you weigh each particle according the
	// scan, instead of multiplying the probabilities, because: log(a*b) = log(a)+log(b).

	// TODO: here comes your code
	calculateDistanceMap(map);
	for (int x = 0; x < likelihoodFieldWidth; x++) {
		for (int y = 0; y < likelihoodFieldHeight; y++) {
			int idx =  x + y * likelihoodFieldWidth;
			double p_hit = Util::gaussian(distMap[idx],sigma * this->likelihoodFieldResolution, 0);
			this->likelihoodField[idx] = log(p_hit *(1-zRand)+ zRand);
		}
	}
	
	ROS_INFO("...DONE creating likelihood field!");
}

void ParticleFilter::calculateDistanceMap(const nav_msgs::OccupancyGrid& map) {
	// calculate distance map = distance to nearest occupied cell
	distMap = new double[likelihoodFieldWidth * likelihoodFieldHeight];
	int occupiedCellProbability = 90;
	// initialize with max distances
	for (int x = 0; x < likelihoodFieldWidth; x++) {
		for (int y = 0; y < likelihoodFieldHeight; y++) {
			distMap[x + y * likelihoodFieldWidth] = 32000.0;
		}
	}
	// set occupied cells next to unoccupied space to zero
	for (int x = 0; x < map.info.width; x++) {
		for (int y = 0; y < map.info.height; y++) {
			if (map.data[x + y * map.info.width] >= occupiedCellProbability) {
				bool border = false;
				for (int i = -1; i <= 1; i++) {
					for (int j = -1; j <= 1; j++) {
						if (!border && x + i >= 0 && y + j >= 0 && x + i
								< likelihoodFieldWidth && y + j
								< likelihoodFieldHeight && (i != 0 || j != 0)) {
							if (map.data[x + i + (y + j) * likelihoodFieldWidth]
									< occupiedCellProbability && map.data[x + i
									+ (y + j) * likelihoodFieldWidth] >= 0)
								border = true;
						}
						if (border)
							distMap[x + i + (y + j) * likelihoodFieldWidth]
									= 0.0;
					}
				}
			}
		}
	}
	// first pass -> SOUTHEAST
	for (int x = 0; x < likelihoodFieldWidth; x++)
		for (int y = 0; y < likelihoodFieldHeight; y++)
			for (int i = -1; i <= 1; i++)
				for (int j = -1; j <= 1; j++)
					if (x + i >= 0 && y + j >= 0 && x + i
							< likelihoodFieldWidth && y + j
							< likelihoodFieldHeight && (i != 0 || j != 0)) {
						double v = distMap[x + i + (y + j)
								* likelihoodFieldWidth] + ((i * j != 0) ? 1.414
								: 1);
						if (v < distMap[x + y * likelihoodFieldWidth]) {
							distMap[x + y * likelihoodFieldWidth] = v;
						}
					}

	// second pass -> NORTHWEST
	for (int x = likelihoodFieldWidth - 1; x >= 0; x--)
		for (int y = likelihoodFieldHeight - 1; y >= 0; y--)
			for (int i = -1; i <= 1; i++)
				for (int j = -1; j <= 1; j++)
					if (x + i >= 0 && y + j >= 0 && x + i
							< likelihoodFieldWidth && y + j
							< likelihoodFieldHeight && (i != 0 || j != 0)) {
						double v = distMap[x + i + (y + j)
								* likelihoodFieldWidth] + ((i * j != 0) ? 1.414
								: 1);
						if (v < distMap[x + y * likelihoodFieldWidth]) {
							distMap[x + y * likelihoodFieldWidth] = v;
						}
					}
}

double* ParticleFilter::getLikelihoodField(int& width, int& height,
		double& resolution) {
	width = this->likelihoodFieldWidth;
	height = this->likelihoodFieldHeight;
	resolution = this->likelihoodFieldResolution;

	return this->likelihoodField;
}

/**
 *  A generic measurement integration method that invokes some specific observation model.
 *  Maybe in the future, we add some other model here.
 */
void ParticleFilter::measurementModel(
		const sensor_msgs::LaserScanConstPtr& laserScan) {
	likelihoodFieldRangeFinderModel(laserScan);
}

/**
 *  Method that implements the endpoint model for range finders.
 *  It uses a precomputed likelihood field to weigh the particles according to the scan and the map.
 */
void ParticleFilter::likelihoodFieldRangeFinderModel(
		const sensor_msgs::LaserScanConstPtr & laserScan) {
		int scan_count = (laserScan->angle_max - laserScan->angle_min)/laserScan->angle_increment;
		scan_count =laser->ranges.size();
		for (int j = 0; j < numberOfParticles; i++) {
			double prob =0;
			for (int i = 0; i<scan_count; i=i+this->laserSkip)
				{
 
					int idx_x = (this->particleSet[i]->x+laserScan->ranges[i]*cos(laserScan->angle_min+laserScan->angle_increment*i))/this->likelihoodFieldResolution;
					int idx_y = (this->particleSet[i]->y+laserScan->ranges[i]*sin(laserScan->angle_min+laserScan->angle_increment*i))/this->likelihoodFieldResolution;
					bool in_range = (laserScan->ranges[i]>laserScan->range_min && laserScan->ranges[i]<laserScan->range_max && );
					int idx = computeMapIndex(likelihoodFieldWidth,likelihoodFieldHeight, idx_x, idx_y);
					bool in_map = idx >=0 && idx<likelihoodFieldHeight * likelihoodFieldWidth;
					if(in_map && in_range)
						{
							this->particleSet[i]->weight += this->likelihoodField[idx];
						}

					else
						{
							// Add a negative weight as punishment
							particleWeight +=log(1E-5);

						}	
				}


}

void ParticleFilter::setMotionModelOdometry(double alpha1, double alpha2,
		double alpha3, double alpha4) {
	this->odomAlpha1 = alpha1;
	this->odomAlpha2 = alpha2;
	this->odomAlpha3 = alpha3;
	this->odomAlpha4 = alpha4;

}

/**
 *  A generic motion integration method that invokes some specific motion model.
 *  Maybe in the future, we add some other model here.
 */
void ParticleFilter::sampleMotionModel(double oldX, double oldY,
		double oldTheta, double newX, double newY, double newTheta) {
	sampleMotionModelOdometry(oldX, oldY, oldTheta, newX, newY, newTheta);

}

/**
 *  Method that implements the odometry-based motion model.
 */
void ParticleFilter::sampleMotionModelOdometry(double oldX, double oldY,
		double oldTheta, double newX, double newY, double newTheta) {
	// TODO: here comes your code
		if(isnan(oldX) || isnan(oldY) || isnan(oldTheta)|| isnan(newTheta) || isnan(newX) || isnan(newY))
		return ;

	//Be careful! the implementaion of motion model is based on all the particles!
	//Knowing what u r doing!
	double delta_trans = std::sqrt((newX-oldX)* (newX-oldX)+(newY-oldY)*(newY-oldY));
	double delta_rot1 = atan2(newY- oldY, newX - oldX);
	double delta_rot2 = (newTheta - oldTheta) - delta_rot1;
	//TODO: should we unwarp the angle
	

 	for (int i = 0; i < numberOfParticles; i++) {
			//Particle* p;
		double delta_rot1_hat = delta_rot1 + Util::gaussianRandom(0, std::abs(this->odomAlpha1*delta_rot1) + this->odomAlpha2 * delta_trans);
		delta_rot1_hat = Util::normalizeTheta(delta_rot1_hat);
		double delta_trans_hat = delta_trans + Util::gaussianRandom(0,this->odomAlpha2 * delta_trans+ this->odomAlpha4*std::abs(delta_rot1 + delta_rot2));
		double delta_rot2 = delta_rot2 + Util::gaussianRandom(0,this->odomAlpha1 * std::abs(delta_rot2)+ this->odomAlpha2 * delta_trans); 
		delta_rot2 = Util::normalizeTheta(delta_rot2);

		this->particleSet[i]->x += delta_trans_hat * cos(this->particleSet[i]->theta + delta_rot1);
		this->particleSet[i]->y += delta_trans_hat * sin(this->particleSet[i]->theta + delta_rot1);
		this->particleSet[i]->theta += delta_rot1_hat + delta_rot2_hat;
		this->particleSet[i]->theta = Util::normalizeTheta(this->particleSet[i]->theta);
		}
}

/**
 *  The stochastic importance resampling.
 */
void ParticleFilter::resample() {

	std::vector<Particle*> newSet;
	double weight_sum = 0;
	int numParticles = getNumberOfParticles();
	double randNumber = Util::uniformRandom(0, 1.0f/(double)particleSet.size() );

	// Cumulative distribution function for resampling.
	double cdf[numParticles];

	cdf[0] = (this->particleSet[0]->weight / this->sumOfParticleWeights);
	// Fill it up with the sum of the normalized weights <= less then i for each i.
	for(int i=1; i < numParticles; i++){
		cdf[i] = cdf[i-1] + this->particleSet[i]->weight / this->sumOfParticleWeights;
	}

	double curThreshold = randNumber;
	for(int j=0, i=0; j < particleSet.size(); j++){

		// Skip until the next threshold reached.
		while(curThreshold > cdf[i]){
			i++;
		}
		// We didn't update the weight as 1/n, since it yields bad results
		Particle* p = new Particle(this->particleSet[i]); 
		newSet.push_back(p);

		// Increment to the next threshold.
		curThreshold += 1.0f/(double)particleSet.size();
	}

	this->particleSet = newSet;
	//TODO: discuss if we need to delete the previous particle

}

Particle* ParticleFilter::getBestHypothesis() {
	double max_weight = 0;
 		for (int i = 0; i < numberOfParticles; i++) {
		if (this->particleSet[i]->weight>max_weight)
		{
			this->bestHypothesis = this->particleSet[i];
			max_weight = this->particleSet[i]->weight;
		}
	}
	return this->bestHypothesis;
}

// added for convenience
int ParticleFilter::computeMapIndex(int width, int height, int x,
		int y) {
	return x + y * width;
}

