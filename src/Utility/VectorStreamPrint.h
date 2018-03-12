#pragma once

template<class T>
ostream & operator<<( ostream & op, std::vector<T> & vec)
{
	for (int i = 0; i < vec.size(); i++)
		op << vec[i] << " ";

	return op;
}

template<class T>
ostream & operator<<( ostream & op, std::vector< std::vector<T> > & vec)
{
	for (int i = 0; i < vec.size(); i++)
		for (int j = 0; j < vec.size(); j++)
			op << vec[i][j] << " ";

	return op;
}

template<class T>
ostream & operator<<( ostream & op, std::vector< std::vector< std::vector<T> > > & vec)
{
	for (int i = 0; i < vec.size(); i++)
		for (int j = 0; j < vec.size(); j++)
			for (int k = 0; k < vec.size(); k++)
				op << vec[i][j][k] << " ";

	return op;
}

template<class T>
ostream & operator<<( ostream & op, std::vector< std::vector< std::vector< std::vector<T> > > > & vec)
{
	for (int i = 0; i < vec.size(); i++)
		for (int j = 0; j < vec.size(); j++)
			for (int k = 0; k < vec.size(); k++)
				for (int l = 0; l < vec.size(); l++)
					op << vec[i][j][k] << " ";

	return op;
}