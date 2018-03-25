#pragma once

#include "../FDV/Node.h"
#include "../Utility/Tensor.h"

#include <fstream>
#include <string>
#include <vector>
// sensors report conservation variable at a given sensor location (currently binds to nearest node) each iteration
// Note: ideally i'd like to use program_options and be consistent with other input, but it looked convoluted to force this kind of syntax

struct Sensor
{
	Sensor(std::string sensorlist, int ndim, std::vector<Node *> & nodeList)
	{
		std::ifstream sensors(sensorlist.c_str());
		std::string name;
		Tensor<double, 1> x(ndim);
		if (ndim == 1)
		{
			std::cout << "1 dimensional sensors are unsupported" << std::endl;
		}
		if (ndim == 2)
		{
			while (sensors >> name >> x(0) >> x(1))
			{
				// process pair (a,b)
				cout << "processing " << name << " " << x << endl;
				sensorList.push_back(nodeList[findSensor(x, nodeList)]);
				streams.emplace_back(std::ofstream{ name + ".sensor" });
			}
		}
		else if (ndim == 3)
		{
			std::cout << "3 dimensional sensors are unsupported" << std::endl;
		}

		// make headers in stream files
		for (int i = 0; i < sensorList.size(); i++)
		{
			streams[i] << "Iteration U V Pressure Temperature" << endl;
		}
	}

	~Sensor()
	{
		for (int i = 0; i < sensorList.size(); i++)
		{
			streams[i].close();
		}
	}

	void update(int j)
	{
		for (int i = 0; i < sensorList.size(); i++)
		{
			streams[i] << j <<  " " << sensorList[i]->v(0) << " " << sensorList[i]->v(1) << " " << sensorList[i]->p << " " << sensorList[i]->T << std::endl;
		}
	}

private:
	int findSensor(Tensor<double, 1> x, std::vector<Node *> & nodeList)
	{
		unsigned int nearest_neighbor = 0;
		double distance = std::numeric_limits<double>::max();

		for (int i = 0; i < nodeList.size(); i++)
		{
			Tensor<double, 1> delta = x - nodeList[i]->x;
//			double d = mag(x - nodeList[i]->x);					// FIXME don't need the sqrt in there. Minor penalty here
			double d = mag(delta);
			if (d < distance)
			{
				nearest_neighbor = i;
				distance = d;
			}
		}
		cout << "nearest neighbor to " << x << " is node " << nearest_neighbor << ", " << nodeList[nearest_neighbor]->x << endl;

		return nearest_neighbor;
	}
	std::vector<std::ofstream> streams;
	std::vector<Node *> sensorList;
};
