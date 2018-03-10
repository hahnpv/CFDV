#pragma once

template<class Elements, class Nodes>
struct AdaptiveRefineList
{
	// temp
	int n;
	AdaptiveRefineList(std::vector<int> & elelist, std::vector<int> & refine_level, std::vector<Elements *> & elements, std::vector<Nodes *> & nodes, int nnod, int nface)
		: elelist(elelist),
		  refine_level(refine_level),
		  elements(elements),
		  nodes(nodes),
		  nnod(nnod),
		  e(elements, nodes, nnod, nface)
	{
		int size = elements.size();
		// int n = 0;
		n = 0;
		for (int i=0; i < size; i++)
		{
			if (elements[i]->refine_level > n)
				n = elements[i]->refine_level;
		}

		refine_level.resize( n+1+1);				// add 1 for initial, add 1 for refine_level of zero
		refine_level[0] = 0;						// start at zero

		FDVParam<Element *> fdvparam;			    /// Calculate the FDV parameters
		double s1, s2, s3, s4;
		for (int i=0; i < size; i++)				// aggregate nodes that need refining. 
		{
			fdvparam(elements[i], s1, s2, s3, s4);	// NOTE used to cache the s values in elements themselves but no longer... was it for RAM or performance?
//			if (elements[i]->s1 > 0.8)
			if (s1 > 0.8)
			{
				refine_list.push_back( i);
			}
		}

		int rsize = refine_list.size();
		for (int i = 0; i < rsize; i++)				// Recursively find dependancies on the original elements to be refined.
		{
			if ( elements[ refine_list[i] ]->refine_level != 0)
			{
				int rl = elements[ refine_list[i] ]->refine_level;
				int el = elements[ refine_list[i] ]->number;
				find_dependancies( rl, el );
			}
		}

		remove_redundant_elements();				// remove redunant elements

		for (int i = 1; i <= n+1; i++)				// now parse list by lvl of refine.
		{
			refine_level[i] = refine_level[i-1];

			for (int j = 0; j < refine_list.size(); j++)
			{
				if ( elements[ refine_list[j] ]->refine_level == i-1 )
				{
					elelist.push_back( refine_list[j] );
					refine_level[i]++;
				}
			}
		}
	};

	void find_dependancies(int refine_level, int elno)
	{
		for (int j = 0; j < e.elements[ elno ].f.size(); j++)
		{
			if ( !e.elements[ elno ].f[j].boundary)
			{
				for( int k = 0; k < e.elements[ elno ].f[j].e.size(); k++)
				{
					if ( e.elements[ elno ].f[j].e[k]->refine_level < refine_level )	// if delta is >1, need multiple passes, but do that later.
					{
						refine_list.push_back( e.elements[ elno ].f[j].e[k]->number );

						int rl = elements[ k ]->refine_level;
						int el = elements[ k ]->number;

						if ( rl > 0)
						{
							find_dependancies( elements[k]->refine_level, elements[ k]->number );
						}
					}
				}
			}
		}

		// add edge dependancies
		for (int j = 0; j < e.elements[ elno ].edge.size(); j++)
		{
			for( int k = 0; k < e.elements[ elno ].edge[j].e.size(); k++)
			{
				if ( e.elements[ elno ].edge[j].e[k]->refine_level < refine_level )	// if delta is >1, need multiple passes, but do that later.
				{
					refine_list.push_back( e.elements[ elno ].edge[j].e[k]->number );	// shouldn't need to walk dependancies
/*																						// but may for very high lvl refines?
					int rl = elements[ k ]->refine_level;
					int el = elements[ k ]->number;

					if ( rl > 0)
					{
						find_dependancies( elements[k]->refine_level, elements[ k]->number );
					}
*/				}
			}
		}

	}

	void remove_redundant_elements()
	{
		for (int i = 0; i < refine_list.size(); i++)
		{
			for (int j = refine_list.size()-1; j > i; j--)
			{
				if (refine_list[i] == refine_list[j])
				{
					refine_list.erase(refine_list.begin() + j, refine_list.begin() + j + 1);
				}
			}
		}
	}

	int nnod;
	std::vector<int> & elelist;					/// Return product, element refinement listing
	std::vector<int> & refine_level;			/// Return product, levels of refinement breakdown of elelist (antiquated eventually?)
	std::vector<Elements *> & elements;			/// Element structure
	std::vector<Nodes *> & nodes;				/// Node structure
	std::vector<int> refine_list;				///	Aggregated nodes that need refining, prior to classification
	ElementAssociation<Element, Node> e;		/// Element association, correlates adjacent elements and nodes
};

