/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h>
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"
#define DIV_VALUE_LIMIT 0.0001 

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[])
{
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of
	//   x, y, theta and their uncertainties from GPS) and all weights to 1.
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	default_random_engine randm;
	
	//number of particles

	num_particles = 120;

	// This line creates a normal (Gaussian) distribution for x.
	normal_distribution<double> distx(x, std[0]);
 
	// TODO: Create normal distributions for y and psi.
	normal_distribution<double> disty(y, std[1]);

	normal_distribution<double> dist_psi(theta, std[2]);

	for (int i = 0; i < num_particles; i++)
	{
		Particle particle;
		particle.x = distx(randm);
		particle.id = i;
		particle.y = disty(randm);
		particle.theta = dist_psi(randm);
		particle.weight = 1;
		particles.push_back(particle);
		weights.push_back(1);
	}
	is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate)
{
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
	default_random_engine randm;

	for (int i = 0; i < num_particles; i++)
	{
		if (fabs(yaw_rate) > DIV_VALUE_LIMIT)
		{
			particles[i].x += velocity / yaw_rate * (sin(particles[i].theta + (yaw_rate * delta_t)) - sin(particles[i].theta));
			particles[i].y += velocity / yaw_rate * (cos(particles[i].theta) - cos(particles[i].theta + (yaw_rate * delta_t)));
			particles[i].theta += yaw_rate * delta_t;
		}
		else
		{
			particles[i].x += velocity * delta_t * cos(particles[i].theta);
			particles[i].y += velocity * delta_t * sin(particles[i].theta);
		}
		normal_distribution<double> distx(particles[i].x, std_pos[0]);

		// TODO: Create normal distributions for y and psi.
		normal_distribution<double> disty(particles[i].y, std_pos[1]);

		normal_distribution<double> dist_psi(particles[i].theta, std_pos[2]);

		particles[i].x = distx(randm);
		particles[i].y = disty(randm);
		particles[i].theta = dist_psi(randm);
	}
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs> &observations)
{
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to
	//   implement this method and use it as a helper during the updateWeights phase.
	for (int i = 0; i < observations.size(); i++)
	{
		double mdist = numeric_limits<double>::max();

		for (int j = 0; j < predicted.size(); j++)
		{
			double _dist = dist(observations[i].x, observations[i].y, predicted[j].x, predicted[j].y);
			if (_dist < mdist)
			{
				mdist = _dist;
				observations[i].id = predicted[j].id;
			}
		}
	}
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
								   std::vector<LandmarkObs> observations, Map map_landmarks)
{
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html
	weights.clear();

	for (int i = 0; i < num_particles; i++)
	{
		vector<LandmarkObs> transformed_observations;

		for (int j = 0; j < observations.size(); j++)
		{
			LandmarkObs transformed_observations_temp;
			transformed_observations_temp.x = observations[j].x * cos(particles[i].theta) - observations[j].y * sin(particles[i].theta) + particles[i].x;
			transformed_observations_temp.y = observations[j].x * sin(particles[i].theta) + observations[j].y * cos(particles[i].theta) + particles[i].y;
			transformed_observations_temp.id = j;
			transformed_observations.push_back(transformed_observations_temp);
		}

		std::vector<LandmarkObs> predicted_landmarks;

		for (int k = 0; k < map_landmarks.landmark_list.size(); k++)
		{
			double dis_temp = dist(particles[i].x, particles[i].y, map_landmarks.landmark_list[k].x_f, map_landmarks.landmark_list[k].y_f);
			if (dis_temp <= sensor_range)
			{
				LandmarkObs landmark;
				landmark.id = map_landmarks.landmark_list[k].id_i;
				landmark.x = map_landmarks.landmark_list[k].x_f;
				landmark.y = map_landmarks.landmark_list[k].y_f;
				predicted_landmarks.push_back(landmark);
			}
		}

		dataAssociation(predicted_landmarks, transformed_observations);

		double weight = 1.0;

		for (int j = 0; j < transformed_observations.size(); j++)
		{
			double meas_x, meas_y, mu_x, mu_y = 0.0;

			meas_x = transformed_observations[j].x;
			meas_y = transformed_observations[j].y;

			for (int k = 0; k < predicted_landmarks.size(); k++)
			{
				if (predicted_landmarks[k].id == transformed_observations[j].id)
				{
					mu_x = predicted_landmarks[k].x;
					mu_y = predicted_landmarks[k].y;
				}
			}

			long double multiplr = exp(-0.5 * (pow(meas_x - mu_x, 2.0) * std_landmark[0] + pow(meas_y - mu_y, 2.0) * std_landmark[1])) / sqrt(2.0 * M_PI * std_landmark[0] * std_landmark[1]);

			if (multiplr > 0)
			{
				weight *= multiplr;
			}
		}
		weights.push_back(weight);
		particles[i].weight = weight;
	}
}

void ParticleFilter::resample()
{
	// TODO: Resample particles with replacement with probability proportional to their weight.
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

	std::discrete_distribution<int> dist(weights.begin(), weights.end());
	default_random_engine randm;
	std::vector<Particle> resampled_particle;
	for (int i = 0; i < num_particles; i++)
	{
		resampled_particle.push_back(particles[dist(randm)]);
	}
	particles = resampled_particle;
}

Particle ParticleFilter::SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y)
{
	//particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
	// associations: The landmark id that goes along with each listed association
	// sense_x: the associations x mapping already converted to world coordinates
	// sense_y: the associations y mapping already converted to world coordinates

	//Clear the previous associations
	particle.associations.clear();
	particle.sense_x.clear();
	particle.sense_y.clear();

	particle.associations = associations;
	particle.sense_x = sense_x;
	particle.sense_y = sense_y;

	return particle;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
	copy(v.begin(), v.end(), ostream_iterator<int>(ss, " "));
	string s = ss.str();
	s = s.substr(0, s.length() - 1); // get rid of the trailing space
	return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
	copy(v.begin(), v.end(), ostream_iterator<float>(ss, " "));
	string s = ss.str();
	s = s.substr(0, s.length() - 1); // get rid of the trailing space
	return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
	copy(v.begin(), v.end(), ostream_iterator<float>(ss, " "));
	string s = ss.str();
	s = s.substr(0, s.length() - 1); // get rid of the trailing space
	return s;
}
