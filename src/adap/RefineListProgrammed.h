#pragma once

void refine_list_programmed(std::vector<int> & elelist, std::vector<int> & refine_level, std::vector<Element *> & elements, int iter)
{
		// try and make a n-level refine
	//	std::vector<int> elelist;			// elements to be refined
	//	std::vector<int> refine_level;		// breakpoints for each level refine
		int size = elements.size();
		int n = 0;
		for (int i=0; i < size; i++)
		{
			if (elements[i]->refine_level > n)
				n = elements[i]->refine_level;
		}

		refine_level.resize( n+1+1);				// add 1 for initial, add 1 for refine_level of zero
		refine_level[0] = 0;						// start at zero

		if ( iter != 4)
		{
		for (int j = 1; j <= n+1; j++)
		{
			refine_level[j] = refine_level[j-1];
			/// TEST
			if ( j == 1 && iter == 4)				// test theory
			{
				// 62, 92, 122 ... < 1320.	
				for (int i = 62; i < 1320; i += 30)
				{	
					elelist.push_back( i);
					refine_level[j]++;
				}
			}
			/// TEST
			for (int i = 0; i < size; i++)			// works so long as we append elements, not insert.
			{
				if (elements[i]->s1 > 0.8 && elements[i]->refine_level == j-1)
				{
					elelist.push_back( i);
					refine_level[j]++;
				}
			}
		}
		}
		else if (iter == 4)
		{
			// tweaked list until we have mapping functionale
			// 0. extend outer bl
			// 62, 92, 122 ... < 1320.	
			// 62, 63, 64...90-1
			for (int i = 62; i < 90; i ++)
			{
				elelist.push_back( i);
				refine_level[1]++;			
			}
			for (int i = 92; i < 1320; i += 30)
			{	
				elelist.push_back( i);
				refine_level[1]++;
			}

			// 1. extend inner bl
			refine_level[2] = refine_level[1];
			for (int i = 0; i < elements.size(); i++)
			{
				if ( elements[i]->refine_level == 1)
				{
					if ( i == 90)
					{
						cout << "hit 90 on the second refine ???" << endl;
						cin.get();
					}
					elelist.push_back( i);
					refine_level[2]++;
				}
			}

			refine_level[3] = refine_level[2];
			// 2. third-level refine on plate
			// all eles with s1>0.8 && refine_level == 1 (couple rows near plate)
			for (int i = 0; i < elements.size(); i ++)
			{	
				if (elements[i]->refine_level == 2)
				{
					elelist.push_back( i );
					refine_level[3]++;
				}
			}
		}
}

