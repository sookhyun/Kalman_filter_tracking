/*!
 *  \file               CellularAutomaton_v1.C
 *  \brief              Implementation of Cellular Automaton for Silicon Tracking detectors
 *  \date    	 	6/4/2018
 *  \author             Sookhyun Lee <dr.sookhyun.lee@gmail.com>
 */


#include "CellularAutomaton_v1.h"
#include "HelixHoughSpace_v1.h"

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/LU>

#include <float.h>
#include <sys/time.h>
#include <cmath>
#include <iostream>

ClassImp(CellularAutomaton_v1)

using namespace std;
using namespace Eigen;

//#define _DEBUG_
// if _FULL_TEST_ is not chosen, only triplets are translated into Track3D
#define _FULL_TEST_


CellularAutomaton_v1::CellularAutomaton_v1(std::vector<Track3D>& input_tracks, std::vector<float>& detector_radii, std::vector<float>& detector_materials)
	:
	_hough_space(NULL),
	_kalman(NULL),
	in_tracks(std::vector<Track3D>()),
	ca_tracks(std::vector<Track3D>()),
	ca_track_states(std::vector<HelixTrackState>()),
	temp_combo(std::vector<unsigned int>()),
	combos(std::set<std::vector<unsigned int> >()),
	layer_sorted(std::vector<std::vector<Cluster3D> >()),
	lbegin(0),
	lend(2),
	nlayers(lend-lbegin+1),
	rlayers(8),
	allowed_missing_inner_hits(0),
	ca_cos_ang_cut(0.985),
	ca_chi2_cut(2.0),
	ca_chi2_layer_cut(2.0),
	ca_phi_cut(M_PI/144.),
	ca_z_cut(20),
	ca_zproj_phi_cut(M_PI/144*1.2),
	ca_zproj_z_cut(1.),
	ca_pproj_phi_cut(M_PI/144*2.),
	ca_pproj_z_cut(0.5),
	ca_dcaxy_cut(0.02),
	fast_chi2_cut_max(FLT_MAX),
	fast_chi2_cut_par0(12.),
	fast_chi2_cut_par1(0.),
	_pt_rescale(1.0),
	_mag_field(1.4),
	_detector_radii(std::vector<float>()),
	_detector_scatter(std::vector<float>()),
	_detector_materials(std::vector<float>()),
	_integrated_scatter(std::vector<float>()),
	_ladder_type(std::vector<int>()),
	_hits_used(std::map<unsigned int, bool>()),
	_hits_map(std::map<unsigned int,Cluster3D>()),
	_layer_clusterid(std::map<int,unsigned int>()),
	CAtime(0.),
	KALtime(0.),
	forward(false),
	remove_hits(false),
	remove_inner_hits(false),
	process_mode(1),
	seeding_mode(false),
	use_projection(false),
	var(999),
	verbosity(10)
{
	set_input_tracks(input_tracks);
	set_detector_radii(detector_radii);
	set_detector_materials(detector_materials);
}

void CellularAutomaton_v1::Reset(){

        if (_hough_space) delete _hough_space;
        if (_kalman) delete _kalman;
        in_tracks.clear();
        ca_tracks.clear();
        ca_track_states.clear();
        
        temp_combo.clear();
	
        _detector_radii.clear();
        _detector_scatter.clear();
        _detector_materials.clear();
        _integrated_scatter.clear();

}


int CellularAutomaton_v1::run(std::vector<Track3D>& output_tracks, std::vector<HelixTrackState>& output_track_states, std::map<unsigned int, bool>& hits_used)
{
	if (remove_hits){
	  	if(verbosity>10){
 	   	cout<< "hits_used size "<< hits_used.size()<<endl;
		cout<< "hits_used size "<< hits_used.size()<<endl;
	  	}
	  _hits_used.swap(hits_used);
	}
	int nintt = (int)_ladder_type.size();
	if (nintt>0) use_projection=true;
	if(verbosity>10) cout<<"process mode "<<process_mode<<" seeding mode "<<seeding_mode<<" projection_mode "<<use_projection<<endl;

	int code = 0;
	cout<< "CellularAutomaton:: initializing..."<<endl;
	code = init();
	cout<<code<<endl;
	if (!code) 
        {
          cout << PHWHERE << "::Error - Initialization failed. " << endl;
          exit(1);
        }
	code = 0;
	cout<<"CellularAutomaton:: processing tracks... "<<endl;
	code = process_tracks();
	if (!code)		
        {
          cout << PHWHERE << "::Error - Processing tracks failed. " << endl;
          exit(1);
        }
	code = 0;
	cout<<"CellularAutomaton:: outputting ca tracks..."<<endl;
	code = get_ca_tracks(output_tracks, output_track_states);
	if (!code)
        {
          cout << PHWHERE << "::Error - Outputting tracks failed. " << endl;
          exit(1);
        }

	if(remove_hits){
	_hits_used.swap(hits_used);		
		if(verbosity>10){
	        cout<< "hits_used size "<< hits_used.size()<<endl;
	        cout<< "_hits_used size "<<_hits_used.size()<<endl;
		}
	}

	for (unsigned int i = 0; i<in_tracks.size(); ++i) in_tracks[i].reset();
	in_tracks.clear();
	for (unsigned int i = 0; i<ca_tracks.size(); ++i) ca_tracks[i].reset();
	ca_tracks.clear();
	ca_track_states.clear();
	layer_sorted.clear();
	return 1;
}

int CellularAutomaton_v1::init()
{
        if (!_hough_space)
        {
          cout << PHWHERE << "::Error - Hough Space is not set. " << endl;
          exit(1);
        }

        if (!_detector_radii.size())
        {
          cout << PHWHERE << "::Error - Detector radii are not set" << endl;
          exit(1);
        }

        if (!_detector_materials.size())
        {
          cout << PHWHERE << "::Error - Detector materials are not set" << endl;
          exit(1);
        }

	nlayers = lend-lbegin+1;
	temp_combo.clear();
	temp_combo.assign(nlayers,0);

	combos.clear();

        std::vector<Cluster3D> one_layer;
        layer_sorted.clear();
	if (process_mode==2){
	layer_sorted.assign(nlayers+3, one_layer);
	} else {
        layer_sorted.assign(nlayers, one_layer);
	}
	if (process_mode==2) sort_clusters_by_layer();
	ca_tracks.clear();
	ca_track_states.clear();
	set_cylinder_kalman();
	
	return 1;
}


void CellularAutomaton_v1::set_hough_space(HelixHoughSpace* hough_space) {

  	_hough_space = hough_space->Clone();
}

void CellularAutomaton_v1::set_mag_field(float mag_field) {
	_mag_field = mag_field;
}

void CellularAutomaton_v1::set_pt_rescale(float pt_rescale){
	_pt_rescale = pt_rescale;
}

void CellularAutomaton_v1::set_detector_radii(std::vector<float>& radii)
{
  for (unsigned int i = 0; i < radii.size(); ++i) {
    _detector_radii.push_back(radii[i]);
  }
}

void CellularAutomaton_v1::set_detector_materials(std::vector<float>& materials)
{
  for (unsigned int i = 0; i < materials.size(); ++i) {
    _detector_scatter.push_back(1.41421356237309515 * 0.0136 *
                               sqrt(3. * materials[i]));
    _detector_materials.push_back(3. * materials[i]);
  }

  _integrated_scatter.assign(_detector_scatter.size(), 0.);
  float total_scatter_2 = 0.;
  for (unsigned int l = 0; l < _detector_scatter.size(); ++l) {
    total_scatter_2 += _detector_scatter[l] * _detector_scatter[l];
    _integrated_scatter[l] = sqrt(total_scatter_2);
  }
}

void CellularAutomaton_v1::set_cylinder_kalman(){
  	_kalman =
      	new HelixKalmanFilter(_detector_radii, _detector_materials, _mag_field);
}

void CellularAutomaton_v1::set_input_tracks(std::vector<Track3D>& input_tracks) 
{
	in_tracks = input_tracks;
	if (verbosity>0) cout<<"Setting input tracks : size = " << in_tracks.size()<<endl; 

}

void CellularAutomaton_v1::sort_clusters_by_layer(){

	// 'extend triplet' mode
        for (std::map<unsigned int,Cluster3D>::iterator jt = _hits_map.begin();
                        jt!= _hits_map.end();
                        ++jt) {

                Cluster3D hit = jt->second;
		unsigned int hitlayer = (unsigned int) hit.get_layer(); 
		if (hitlayer>2) hitlayer = (hitlayer-3)/2 +3;
		unsigned int layer = hitlayer - lbegin +3;
                if(verbosity>10) cout<<"layer "<<layer<< endl;
                if (!forward) layer = nlayers+3-layer-1;
                if (layer > (nlayers+3-1)) continue;
                unsigned int min = (layer - allowed_missing_inner_hits);
                if (allowed_missing_inner_hits > layer) {
                min = 0;
                }
                for (unsigned int l = min; l <= layer; ++l) {
                if(l>=0 && l<3) continue; // require first 3 hits, no missing hits allowed 
                layer_sorted[l].push_back(hit);
                if(verbosity>10)cout<<"adding hit in layer "<<l<<endl;
                }
        }

}

int CellularAutomaton_v1::get_ca_tracks(std::vector<Track3D>& output_tracks, std::vector<HelixTrackState>& output_track_states)
{
	// push back new ca processed tracks into _tracks
	   
	if (ca_tracks.size() != ca_track_states.size()) 
	return 0;

	for (unsigned int i = 0; i <ca_tracks.size(); ++i)
	{
	output_tracks.push_back(ca_tracks[i]);
	output_track_states.push_back(ca_track_states[i]);
	}

	if (verbosity>0)cout<<"newly added ca tracks : "<< ca_tracks.size() <<" tracks. "<<endl;
	return 1;
}


int CellularAutomaton_v1::process_tracks()
{

	for (unsigned int i = 0; i < in_tracks.size(); ++i) 
	{ // loop over input tracks
		if (verbosity>1)cout<<"track candidate "<<i<<endl;

		switch(process_mode){

		case 0: // Tracks from Hough  Transform
		process_single_houghbin(in_tracks[i]);
		break;

		case 1: // Construct continuous tracklets (triplet,quadraplet and etc.), forward mode only for now
		process_single_tracklet(in_tracks[i]);
		break;// case 1

		case 2: // Start from a built triplet and extend it with an option of using projected z 
		extend_single_triplet(in_tracks[i]);
		break;

		default:
		return 0;
		}
	}

	return 1;

}


int CellularAutomaton_v1::process_single_tracklet(Track3D& track){ 

        unsigned int cur_seg_size = 0;
        unsigned int next_seg_size = 0;

	std::map<unsigned int, TrackSegment> cur_seg;
	std::map<unsigned int, TrackSegment> next_seg;


//	std::vector<TrackSegment> complete_segments;
//	unsigned int comp_seg_size = 0;

/***
        Convert cluster triplets to a segment with associated kappa & dzdl. 
	Sort out clusters by layers first and allow for non-zero missing layers.
 	l=0 : layer lbegin, l=nlayers-1 : layer lend
***/
	for (unsigned int l = 0; l < nlayers; ++l) { 
		layer_sorted[l].clear();
	}

  	for (unsigned int i = 0; i < track.hits.size(); ++i) {
    		Cluster3D hit = track.hits[i];
                unsigned int hitlayer = (unsigned int) hit.get_layer();
                if (hitlayer>2) hitlayer = (hitlayer-3)/2 +3;
                unsigned int layer = hitlayer;

		//cout<<"layer "<<layer<< endl;
       	 	if (!forward) layer = nlayers-layer-1;
                if (layer > (nlayers-1)) continue;
    		unsigned int min = (layer - allowed_missing_inner_hits);
    		if (allowed_missing_inner_hits > layer) {
      		min = 0;
    		}
    		for (unsigned int l = min; l <= layer; ++l) {
      		layer_sorted[l].push_back(hit);
		//cout<<"adding hit in layer "<<l<<endl;
		}
	}

	for (unsigned int l = 0; l< 3; ++l){
	if (verbosity>1)cout<<"layer_sorted["<<l<<"].size = "<< layer_sorted[l].size()<<endl;
        	if (layer_sorted[l].size() == 0) {
        	return 0;
    		}
  	}

#ifdef _FULL_TEST_
	float ca_cos_ang_cut_diff = 1. - ca_cos_ang_cut;
	float ca_cos_ang_cut_diff_inv = 1. / ca_cos_ang_cut_diff;
#endif
	float ca_sin_ang_cut = sqrt(1. - ca_cos_ang_cut * ca_cos_ang_cut);

	std::vector<float> inv_layer;
	inv_layer.assign(nlayers, 1.);
	for (unsigned int l = 3; l < nlayers; ++l) {
	inv_layer[l] = 1. / (((float)l) - 2.);
	}

	// l = 3, 4,   5,   6,   7
	//     1, 1/2, 1/3, 1/4, 1/5
	//     1, 1/3, 1/5, 1/7, 1/9         

	float x1, x2, x3;
	float y1, y2, y3;
	float z1, z2, z3;
	float dx1, dx2, dx3;
	float dy1, dy2, dy3;
	float dz1, dz2, dz3;

 	float kappa;
	float dkappa;

	float ux_mid;
 	float uy_mid;
	float ux_end;
	float uy_end;

	float dzdl_1;
	float dzdl_2;
	float ddzdl_1;
	float ddzdl_2;

#ifdef _FULL_TEST_
	float cur_kappa;
	float cur_dkappa;
	float cur_ux;
	float cur_uy;
	float cur_chi2;
	float chi2;
#endif
	unsigned int hit1;
	unsigned int hit2;
	unsigned int hit3;

	TrackSegment temp_segment;
	temp_segment.hits.assign(nlayers, 0);
  	for (unsigned int i = 0; i < layer_sorted[0].size(); ++i) {
		for (unsigned int j = 0; j < layer_sorted[1].size(); ++j) {
			for (unsigned int k = 0; k < layer_sorted[2].size() ; ++k) {

        		x1 = layer_sorted[0][i].get_x();
        		y1 = layer_sorted[0][i].get_y();
       		 	z1 = layer_sorted[0][i].get_z();

        		dx1 = 0.5*sqrt(12.0)*sqrt(layer_sorted[0][i].get_size(0,0));
        		dy1 = 0.5*sqrt(12.0)*sqrt(layer_sorted[0][i].get_size(1,1));
        		dz1 = 0.5*sqrt(12.0)*sqrt(layer_sorted[0][i].get_size(2,2));

        		x2 = layer_sorted[1][j].get_x();
        		y2 = layer_sorted[1][j].get_y();
        		z2 = layer_sorted[1][j].get_z();

        		dx2 = 0.5*sqrt(12.0)*sqrt(layer_sorted[1][j].get_size(0,0));
        		dy2 = 0.5*sqrt(12.0)*sqrt(layer_sorted[1][j].get_size(1,1));
        		dz2 = 0.5*sqrt(12.0)*sqrt(layer_sorted[1][j].get_size(2,2));

        		x3 = layer_sorted[2][k].get_x();
        		y3 = layer_sorted[2][k].get_y();
        		z3 = layer_sorted[2][k].get_z();
        		dx3 = 0.5*sqrt(12.0)*sqrt(layer_sorted[2][k].get_size(0,0));
        		dy3 = 0.5*sqrt(12.0)*sqrt(layer_sorted[2][k].get_size(1,1));
        		dz3 = 0.5*sqrt(12.0)*sqrt(layer_sorted[2][k].get_size(2,2));

        		hit1 = i;
       	 		hit2 = j;
        		hit3 = k;

        		calculate_kappa_tangents(x1, y1, z1, x2, y2, z2, x3, y3,
                                z3, dx1, dy1, dz1, dx2, dy2, dz2,
                                dx3, dy3, dz3, kappa, dkappa,
                                ux_mid, uy_mid, ux_end, uy_end,
                                dzdl_1, dzdl_2, ddzdl_1, ddzdl_2);

#ifdef _DEBUG_
        cout<<"Triplet : "<<endl;
        cout<<"kappa "<<kappa<< " dkappa "<<dkappa<<" ux_mid "<<ux_mid<<" uy_mid "<<" ux_end "<<ux_end<<
        " uy_end " <<uy_end<<" dzdl_1 "<<dzdl_1<<" dzdl_2 "<<dzdl_2<<" ddzdl_1 "<<ddzdl_1<<" ddzdl_2 "<<ddzdl_2<<endl;
#endif

          		temp_segment.chi2 = pow(
              		(dzdl_1 - dzdl_2) /
              		(ddzdl_1 + ddzdl_2 + fabs(dzdl_1 * ca_sin_ang_cut)),2);

         		if (temp_segment.chi2 > ca_chi2_layer_cut) continue;
          		temp_segment.ux = ux_end;
          		temp_segment.uy = uy_end;
          		temp_segment.kappa = kappa;
			//if (temp_segment.kappa > _hough_space->get_kappa_max()) continue;
                      	temp_segment.dkappa = dkappa;

          		temp_segment.hits[0] = hit1;
          		temp_segment.hits[1] = hit2;
          		temp_segment.hits[2] = hit3;
          		temp_segment.n_hits = 3;

          		unsigned int outer_layer =
                	layer_sorted[2][temp_segment.hits[2]].get_layer();
          		if (!forward) outer_layer = nlayers-outer_layer-1;
          		// make sure we have required number of layers with hits
                       	if ((outer_layer - 2) > allowed_missing_inner_hits) continue;
          		// finish up if required number of layers is reached,
                       	if ((nlayers - 3) <= allowed_missing_inner_hits) {
//                	complete_segments.push_back(temp_segment);
//			++comp_seg_size;
          		}


          		if (next_seg.size() == next_seg_size) { // first new segment
                        next_seg.insert(make_pair(next_seg_size, temp_segment));
                        next_seg_size += 1;
          		} else { // next new segments
               			(next_seg)[next_seg_size] = temp_segment;
                		next_seg_size += 1;
          		}

			}
		}
	}

	cur_seg.swap(next_seg);
  	swap(cur_seg_size, next_seg_size);


  	if (verbosity>2) cout<<"number of complete segments : " << comp_seg_size<<endl;
  	if (verbosity>2) cout<<"number of current segments : "<< cur_seg_size<<endl;
/*
  	// copy complete segments over to current segments
       	for (unsigned int i = 0; i < comp_seg_size; ++i) {
//    		if (cur_seg->size() == cur_seg_size) {
//      		cur_seg->push_back(complete_segments[i]);
//      		++cur_seg_size;
//    		} else {
     	 	(*cur_seg)[cur_seg_size] = complete_segments[i];
      		++cur_seg_size;
//    		}
  	}
	cout<<"here 3"<<endl;
  	std::set<unsigned int> comp1;
  	std::set<unsigned int> comp2;

//  	cout<<"number of segments generated "<< cur_seg_size<<endl;
  	if (cur_seg_size==0  || cur_seg_size>10000) return 1;
  	for (unsigned int i = cur_seg_size-1; i>0; --i)
	{
    		if ((*cur_seg)[i].n_hits==0) continue;
    		comp1.clear();
    		temp_combo.assign((*cur_seg)[i].n_hits, 0);
    		for (unsigned int l = 0; l < (*cur_seg)[i].n_hits; ++l) {
    		temp_combo[l] = layer_sorted[l][(*cur_seg)[i].hits[l]].get_id();
    		comp1.insert(temp_combo[l]);
    		}

    		sort(temp_combo.begin(), temp_combo.end());
    		set<vector<unsigned int> >::iterator it = combos.find(temp_combo);
    		if (it != combos.end()) {
     	 	(*cur_seg)[i].n_hits = 0.;
    		}
    		if (combos.size() > 10000) {
      		combos.clear();
    		}
    		combos.insert(temp_combo);

    		for (unsigned int j = i-1; j>=0; --j)
		{
      			comp2.clear();
      			for (unsigned int m = 0; m < (*cur_seg)[j].n_hits; ++m){
      				comp2.insert(layer_sorted[m][(*cur_seg)[j].hits[m]].get_id());
      			};
      			for (std::set<unsigned int>::iterator it=comp1.begin(); it!=comp1.end(); ++it){
      			auto it2 = comp2.find(*it);
      			if (it2 != comp2.end()) comp2.erase(*it2);
      			}
      			if (comp2.empty()) {
        		(*cur_seg)[j].n_hits = 0;
      			}
      			if (j==0) break;
    		}
  	}

	comp1.clear();
	comp2.clear();
*/
  	unsigned int nsegs = cur_seg_size;
  	for (unsigned int i = 0; i<cur_seg_size; i++) {
  		if ((cur_seg)[i].n_hits ==0) --nsegs;
  	}
  	if(verbosity>0) cout<<"number of segments from triplets to be processed "<< nsegs<<endl;


#ifdef _FULL_TEST_

// start over counting after triplets
// To Do: deal with allowing hits missing in main module
//	unsigned int allowed_missing = nlayers - rlayers;
////	cout<<"allowed missing "<< allowed_missing<<endl;

	std::map<unsigned int, unsigned int> missing_layers_map; // segment_id, missing_layers
	std::map<unsigned int, unsigned int> missing_layers_map_next;
	missing_layers_map.clear();
	missing_layers_map_next.clear();
	for (unsigned int n = 0 ; n < cur_seg_size; ++n) missing_layers_map.insert(make_pair(n,0));
	unsigned int added_next_segments = 0;

 	// Extend segments, start from 4th layer, does not have to be INTT if we are doing TPC
 	for (unsigned int l=3; l < nlayers; ++l ) {
		//cout<<"nlayers "<<nlayers<<" l "<<l<<endl;
		next_seg_size = 0;
		next_seg.clear();
		// Loop over current segments
		for (unsigned int iseg = 0; iseg < cur_seg_size; ++iseg) {
			if ((cur_seg)[iseg].n_hits ==0) continue;
			added_next_segments = 0;
			auto search = missing_layers_map.find(iseg);
                        unsigned int missing_layers = search->second;
			//cout<<"segment "<< i <<", missing layer "<<missing_layers<<endl;
			// loop over all hits on layer n
			// if multiple legitimate hits in next layer, duplicate current segment, copy over missing layer and assign new segment_id
			//// ** drop layer sorted keep good hits in layer_sorted, hits[l] : id 
			// keep cluster ids in hits[l] for good hits of a segment
			// for a layer without a good hit, store with hits[l] = 999

			if ( (l-2) < 3){
		        x1 = layer_sorted[l - 2][(cur_seg)[iseg].hits[l - 2]].get_x();
		        y1 = layer_sorted[l - 2][(cur_seg)[iseg].hits[l - 2]].get_y();
        		z1 = layer_sorted[l - 2][(cur_seg)[iseg].hits[l - 2]].get_z();
		        dx1 = 0.5*sqrt(12.0)*sqrt(layer_sorted[l - 2][(cur_seg)[iseg].hits[l - 2]].get_size(0,0));
                        dy1 = 0.5*sqrt(12.0)*sqrt(layer_sorted[l - 2][(cur_seg)[iseg].hits[l - 2]].get_size(1,1));
                        dz1 = 0.5*sqrt(12.0)*sqrt(layer_sorted[l - 2][(cur_seg)[iseg].hits[l - 2]].get_size(2,2));
			}else {
			// get positions from hits_map
			auto searchl2 = _hits_map.find((cur_seg)[iseg].hits[l-2]);
			Cluster3D clusterl2 = searchl2->second;
			x1 = clusterl2.get_x();
			y1 = clusterl2.get_y();
			z1 = clusterl2.get_z();
			dx1 = 0.5*sqrt(12.0)*sqrt(clusterl2.get_size(0,0));
			dy1 = 0.5*sqrt(12.0)*sqrt(clusterl2.get_size(1,1));
			dz1 = 0.5*sqrt(12.0)*sqrt(clusterl2.get_size(2,2));
			}

			if  ( (l-1) < 3) {
        		x2 = layer_sorted[l - 1][(cur_seg)[iseg].hits[l - 1]].get_x();
        		y2 = layer_sorted[l - 1][(cur_seg)[iseg].hits[l - 1]].get_y();
        		z2 = layer_sorted[l - 1][(cur_seg)[iseg].hits[l - 1]].get_z();
                        dx2 = 0.5*sqrt(12.0)*sqrt(layer_sorted[l - 1][(cur_seg)[iseg].hits[l - 1]].get_size(0,0));
                        dy2 = 0.5*sqrt(12.0)*sqrt(layer_sorted[l - 1][(cur_seg)[iseg].hits[l - 1]].get_size(1,1));
                        dz2 = 0.5*sqrt(12.0)*sqrt(layer_sorted[l - 1][(cur_seg)[iseg].hits[l - 1]].get_size(2,2));
			} else {	
			// get positions from hits_map
			auto searchl1 = _hits_map.find((cur_seg)[iseg].hits[l-1]);
			Cluster3D clusterl1 = searchl1->second;
			x2 = clusterl1.get_x();
			y2 = clusterl1.get_y();
			z2 = clusterl1.get_z(); 
                        dx2 = 0.5*sqrt(12.0)*sqrt(clusterl1.get_size(0,0));
                        dy2 = 0.5*sqrt(12.0)*sqrt(clusterl1.get_size(1,1));
                        dz2 = 0.5*sqrt(12.0)*sqrt(clusterl1.get_size(2,2));
			}

                        // convert segment to track to process it through kalman filter or just call fit_track

			Track3D init_track;
			bool fit_layer = (l >= 6);
			if (fit_layer)
			{       //init_track.hits.assign(l, Cluster3D());
                                init_track.hits.assign((cur_seg)[iseg].n_hits, Cluster3D());

                		for (unsigned int ll = 0; ll < (cur_seg)[iseg].n_hits; ++ll) 
				{
                			if (ll<3)
					{
                			init_track.hits[ll] = layer_sorted[ll][(cur_seg)[iseg].hits[ll]];
                			}else 
					{
                			auto search = _hits_map.find((cur_seg)[iseg].hits[ll]);
                			Cluster3D cluster = search->second;
                			init_track.hits[ll] = cluster;
                			}
		                }
			}
			// Now we are on the new layer, layer l !
                       	for (std::map<unsigned int,Cluster3D>::iterator jt = _hits_map.begin();
                        		jt!= _hits_map.end();
                                        ++jt) 
			{
                        	Cluster3D hit3d = jt->second;
                        	unsigned int layer = hit3d.get_layer();
                        	if (layer != l) continue;
//				which_seg = i;
                                hit1 = hit3d.get_id();
				
                                x3 = hit3d.get_x();
                                y3 = hit3d.get_y();
                                z3 = hit3d.get_z();

                                float phi_prev =shift_phi_range(atan2(y2,x2));
                                float phi_cur = shift_phi_range(atan2(y3,x3));
				float phi_diff = phi_cur-phi_prev;

				if (phi_cur< M_PI/2. && phi_prev > 3*M_PI/2.) phi_diff += 2.*M_PI;
				else if (phi_cur>3*M_PI/2 && phi_prev<M_PI/2.) phi_diff -= 2.*M_PI;

				if (!seeding_mode){
                                	if ((fabs(phi_diff)> ca_phi_cut || abs(z3-z2)> ca_z_cut)) continue;
				} else {
					if ((fabs(phi_diff)> ca_phi_cut || abs(z3-z2)> ca_z_cut) && cur_seg_size!=1) continue;
				}
					
//                		auto search = _hits_used.find(hit1);
//                		if(search != _hits_used.end() && search->second ) continue;

				switch (fit_layer){
	
				case false :

        	                dx3 = 0.5*sqrt(12.0)*sqrt(hit3d.get_size(0,0));
                	        dy3 = 0.5*sqrt(12.0)*sqrt(hit3d.get_size(1,1));
                        	dz3 = 0.5*sqrt(12.0)*sqrt(hit3d.get_size(2,2));

       		 		cur_kappa = (cur_seg)[iseg].kappa;
        			cur_dkappa = (cur_seg)[iseg].dkappa;
        			cur_ux = (cur_seg)[iseg].ux;
        			cur_uy = (cur_seg)[iseg].uy;
        			cur_chi2 = (cur_seg)[iseg].chi2;

				chi2=9999;
	        		calculate_kappa_tangents(
        	      			x1, y1, z1, x2, y2, z2, x3, y3, z3,
              				dx1, dy1, dz1, dx2, dy2, dz2, dx3, dy3, dz3,
              				kappa, dkappa, ux_mid, uy_mid, ux_end, uy_end,
              				dzdl_1, dzdl_2, ddzdl_1, ddzdl_2,
              				ca_sin_ang_cut, ca_cos_ang_cut_diff_inv,
              				cur_kappa, cur_dkappa, cur_ux, cur_uy, cur_chi2, chi2);
#ifdef _DEBUG_
        			cout<<"Extended layers for segment "<<which_seg<<endl;
        			cout<<"kappa "<<kappa<< " dkappa "<<dkappa
        			<<" ux_mid "<<ux_mid<<" uy_mid "<<" ux_end "<<ux_end<<" uy_end " <<uy_end
        			<<" dzdl_1 "<<dzdl_1<<" dzdl_2 "<<dzdl_2<<" ddzdl_1 "<<ddzdl_1<<" ddzdl_2 "<<ddzdl_2
        			<<" cur_chi2 "<<cur_chi2<<" chi2 "<<chi2<<" chi2*inv_layer "<<l <<" "<<chi2*inv_layer[l]
        			<<" chi2_cut " <<ca_chi2_layer_cut<<endl;
#endif


			// if segment with good chi2, store it in next segment
        			if (chi2 * inv_layer[l] < ca_chi2_layer_cut ) {
              				temp_segment.chi2 = chi2;
              				temp_segment.ux = ux_end;
              				temp_segment.uy = uy_end;
              				temp_segment.kappa = kappa;
					if (seeding_mode){ // this is a temporaty fix to get high efficiency triplet seeds 
              					if (temp_segment.kappa > _hough_space->get_kappa_max() &&cur_seg_size!=1) continue;
					} else {
						if (temp_segment.kappa > _hough_space->get_kappa_max()) continue;
					} 
					
                        		temp_segment.dkappa = dkappa;
              				for (unsigned int ll = 0; ll < l; ++ll) {
                			temp_segment.hits[ll] = (cur_seg)[iseg].hits[ll];
              				}
              				temp_segment.hits[l] = hit1;
              				//unsigned int outer_layer =
                  			//	layer_sorted[l][temp_segment.hits[l]].get_layer();
              				//if (!forward) outer_layer = nlayers-outer_layer-1;
              				temp_segment.n_hits = l + 1;
              				// finish up if required number of layers is reached
                        		//if ((nlayers - (l + 1)) <= allowed_missing) {
                			//	complete_segments.push_back(temp_segment);
              				//}
              				// make sure we have required number of layers with hits
                        		//if ((outer_layer - l) > allowed_missing) {
                			//continue;
              				//}

		                	// Discard segment if too many missing layers                 		
              				if (next_seg.size() == next_seg_size) { // first new segment
                				next_seg.insert(make_pair(next_seg_size,temp_segment));
                                        	missing_layers_map_next.insert(make_pair(next_seg.size()-1, missing_layers));
                				next_seg_size += 1;
              				} else { // next new segments
					#ifdef _DEBUG_
					cout<<"Next segment size inconsistent "<<endl;
					#endif
                				(next_seg)[next_seg_size] = temp_segment;
				        	missing_layers_map_next.insert(make_pair(next_seg.size()-1, missing_layers));
                				next_seg_size += 1;
              				}
					++added_next_segments;
#ifdef _DEBUG_					
					cout<<"segment "<< iseg << " added segment "<<added_next_segments<<endl;
#endif
         			} // chi2 cut from segment building method -> change to switch - case block
			
				break;

				case true :				

				// fit init_track to get kappa to compare with new kappa
				// copy init_track over to temp_track and add a hit in kalman filter
				Track3D temp_track;
			        temp_track.hits.assign(init_track.hits.size()+1, Cluster3D());

      				for (unsigned int ll = 0; ll < init_track.hits.size(); ++ll) {
        			temp_track.hits[ll] = init_track.hits[ll];
      				}
				temp_track.hits[init_track.hits.size()] = hit3d;

				// track fitting instead of computing from triplets

 				float temp_chi2 = temp_track.fit_track();
//				cout<<"chi2 from fit_track "<<init_chi2 <<endl;			
				if (temp_chi2 != temp_chi2) continue;
		                if (temp_track.kappa != temp_track.kappa )  continue;
                		if (temp_track.z0 != temp_track.z0) continue;

		                HelixTrackState state;
		                state.phi = temp_track.phi;
		                if (state.phi < 0.) {
		                state.phi += 2. * M_PI;
		                }
		                state.d = temp_track.d;
		                state.kappa = temp_track.kappa;
		                state.nu = sqrt(state.kappa);
		                state.z0 = temp_track.z0;
		                state.dzdl = temp_track.dzdl;
		                state.C = Matrix<float, 5, 5>::Zero(5, 5);
		                state.C(0, 0) = pow(0.01, 2.);
		                state.C(1, 1) = pow(0.01, 2.);
		                state.C(2, 2) = pow(0.01 * state.nu, 2.);
		                state.C(3, 3) = pow(0.05, 2.);
		                state.C(4, 4) = pow(0.05, 2.);
		                state.chi2 = 0.;
		                state.position = 0;
		                state.x_int = 0.;
		                state.y_int = 0.;
		                state.z_int = 0.;
			
				// place holder for kalman filter
/*
	                	unsigned int nfits = 0;
                		for (unsigned int h = 0; h < temp_track.hits.size(); ++h) {
                		_kalman->addHit(temp_track.hits[h], state);
                		nfits += 1;
                		cout<<"nfits "<<nfits<<endl;
                		}
                		cout<<"z0 after kalman "<<state.z0<<endl;

                		// fudge factor for non-gaussian hit sizes
                                state.C *= 3.;
                		state.chi2 *= 6.;
                		// kappa cut here drives both efficiency and ghost track rates down at the same time
                                if (!(temp_track.kappa == temp_track.kappa) ) {
                		continue;
                		}
                		if (!(state.chi2 == state.chi2)) {
                		continue;
                		}

*/

              			if (state.chi2 / (2. * ((float)(temp_track.hits.size())) - 5.) < ca_chi2_cut) {
					// translate temp_track into temp_segment (only hit info is saved) and save in next segments
					for (unsigned int ll = 0; ll < l; ++ll) {
                                        temp_segment.hits[ll] = (cur_seg)[iseg].hits[ll];
                                        }
                                        temp_segment.hits[l] = hit1;
					temp_segment.n_hits = l + 1;

                                        if (next_seg.size() == next_seg_size) { // first new segment
                                                next_seg.insert(make_pair(next_seg_size,temp_segment));
                                                missing_layers_map_next.insert(make_pair(next_seg.size()-1, missing_layers));
                                                next_seg_size += 1;
                                        } else { // next new segments
                                        #ifdef _DEBUG_
                                        cout<<"Next segment size inconsistent "<<endl;
                                        #endif
                                                (next_seg)[next_seg_size] = temp_segment;
                                                missing_layers_map_next.insert(make_pair(next_seg.size()-1, missing_layers));
                                                next_seg_size += 1;
                                        }
                                        ++added_next_segments;
#ifdef _DEBUG_
                                        cout<<"segment "<< which_seg << " added segment "<<added_next_segments<<endl;
#endif
                		
                		}// chi2 from track fitting and/or kalman filter : looser cuts than computing kappa from triplets
	
				break;

				}// l>=6 switch	


			}// clusters on layer l		
			
/*
			if (added_next_segments==0){
			 ++missing_layers;
			if (missing_layers <= allowed_missing){
                                if (next_seg->size() == next_seg_size) { // first new segment
                                        next_seg->push_back((*cur_seg)[i]);
					missing_layers_map_next.insert(make_pair(next_seg->size()-1, missing_layers));
                                        next_seg_size += 1;
                                } else { // next new segments
                                        (*next_seg)[next_seg_size] = (*cur_seg)[i];
					missing_layers_map_next.insert(make_pair(next_seg->size()-1, missing_layers));
                                        next_seg_size += 1;
                                }
			}
			}
*/


        	} // current segments
	    	cur_seg.swap(next_seg);
    		swap(cur_seg_size, next_seg_size);
		missing_layers_map.swap(missing_layers_map_next);
		// Current segments now hold extended segments up to layer l, and missing layers info is in missing_layers_map 

	}// next segment


#endif // FULL_TEST

	Track3D temp_track;
	temp_track.hits.assign(nlayers, Cluster3D());

	std::vector<Track3D> best_track;
	std::vector<HelixTrackState> best_track_state;
	float best_chi2 = 9999;
	for (unsigned int i = 0; i< cur_seg_size; ++i) {

		if ((cur_seg)[i].n_hits ==0) continue;

#ifdef _DEBUG_
		cout<<"segment " <<i <<endl;
#endif
      		temp_track.hits.assign((cur_seg)[i].n_hits, Cluster3D());
      		for (unsigned int l = 0; l < (cur_seg)[i].n_hits; ++l) {
			if (l<3){
	        	temp_track.hits[l] = layer_sorted[l][(cur_seg)[i].hits[l]];
			}else {
			auto search = _hits_map.find((cur_seg)[i].hits[l]);
			Cluster3D cluster = search->second;
			temp_track.hits[l] = cluster;
		}
      		}

      		float init_chi2 = temp_track.fit_track();
                if (verbosity>1) cout<<"chi2 from fit_track "<<init_chi2 <<" kappa "<< temp_track.kappa <<endl;

#ifdef _DEBUG_
          	cout	<<" kappa " <<temp_track.kappa <<" phi "<<temp_track.phi<<" d "<<temp_track.d
          		<<" z0 "<<temp_track.z0<<" dzdl "<<temp_track.dzdl<< endl;
#endif
		
		if (seeding_mode){
			if (temp_track.kappa != temp_track.kappa && cur_seg_size!=1) continue;
			if (temp_track.z0 != temp_track.z0 && cur_seg_size != 1) continue;
		} else {
		        if (temp_track.kappa != temp_track.kappa) continue;
                        if (temp_track.z0 != temp_track.z0) continue;
		}

    		HelixTrackState state;
    		state.phi = temp_track.phi;
    		if (state.phi < 0.) {
      		state.phi += 2. * M_PI;
    		}
    		state.d = temp_track.d;
   		state.kappa = temp_track.kappa;
    		state.nu = sqrt(state.kappa);
    		state.z0 = temp_track.z0;
    		state.dzdl = temp_track.dzdl;
    		state.C = Matrix<float, 5, 5>::Zero(5, 5);
    		state.C(0, 0) = pow(0.01, 2.);
    		state.C(1, 1) = pow(0.01, 2.);
    		state.C(2, 2) = pow(0.01 * state.nu, 2.);
    		state.C(3, 3) = pow(0.05, 2.);
    		state.C(4, 4) = pow(0.05, 2.);
    		state.chi2 = 0.;
    		state.position = 0;
    		state.x_int = 0.;
    		state.y_int = 0.;
    		state.z_int = 0.;

    		unsigned int nfits = 0;
    		for (unsigned int h = 0; h < temp_track.hits.size(); ++h) {
      		_kalman->addHit(temp_track.hits[h], state);
      		nfits += 1;
#ifdef _DEBUG_
		cout<<"nfits "<<nfits<<endl;
#endif
    		}

    		if (verbosity>1) cout<<"z0 after kalman "<<state.z0<<endl;
   		// fudge factor for non-gaussian hit sizes
           	state.C *= 3.;
    		state.chi2 *= 6.;
		// kappa cut *here* drives both efficiency and ghost track rates down at the same time
		if (seeding_mode){
    		if (!(temp_track.kappa == temp_track.kappa) && (cur_seg_size !=1)) continue;
    		if (!(state.chi2 == state.chi2) && (cur_seg_size !=1)) continue;
		} else {
                if (!(temp_track.kappa == temp_track.kappa)) continue;
                if (!(state.chi2 == state.chi2)) continue;
		}

/*
      		if (fabs(temp_track.d) > ca_dcaxy_cut) continue;
      		if (fabs(temp_track.z0) > dca_cut) continue;
*/


		// no chi2 cut for tests on triplets
		if (verbosity>1) cout<<"state.chi2 from kalman "<<state.chi2<<endl;
		if (seeding_mode){
    			if (state.chi2 / (2. * ((float)(temp_track.hits.size())) - 5.) > ca_chi2_cut && cur_seg_size !=1) continue;
		} else {
			if (state.chi2 / (2. * ((float)(temp_track.hits.size())) - 5.) > ca_chi2_cut) continue;
		}

    		if (best_chi2 > state.chi2 || (seeding_mode && cur_seg_size==1)){
      			if (!best_track.empty()){
      			best_track.pop_back();
      			best_track_state.pop_back();
      			}
      			best_track.push_back(temp_track);
     			best_track_state.push_back(state);
    		}
  	}

	if (best_track.empty()) return 1;

	if(verbosity >2) cout<<"best_track.size "<<best_track.size()<<endl;
  	ca_tracks.push_back(best_track.back());
  	ca_track_states.push_back(best_track_state.back());
  	if (verbosity>1) cout <<"ca track added, chi2 =  "<< (best_track_state.back().chi2)/(2. * ((float)(temp_track.hits.size())) - 5.) <<" z0 = "<<best_track_state.back().z0<<endl;


	//    if ((remove_hits == true) && (state.chi2 < chi2_removal_cut) &&
      	//        (temp_track.hits.size() >= n_removal_hits)) {
        temp_track = best_track.back();

        if  (remove_hits){
        for (unsigned int i = 0; i < temp_track.hits.size(); ++i) {
                if (!remove_inner_hits && temp_track.hits[i].get_layer()<3) continue;
                auto search = _hits_used.find(temp_track.hits[i].get_id());
                if(search != _hits_used.end())
                {
                _hits_used.find(temp_track.hits[i].get_id())->second = true;
                }
        }
        }

	cur_seg.clear();
	next_seg.clear();
	return 1;
}

int CellularAutomaton_v1::process_single_houghbin(Track3D& track)
{

        std::vector<TrackSegment> segments1;
        std::vector<TrackSegment> segments2;


  std::vector<TrackSegment>* cur_seg = &segments1;
  std::vector<TrackSegment>* next_seg = &segments2;
  unsigned int cur_seg_size = 0;
  unsigned int next_seg_size = 0;

  std::vector<TrackSegment> complete_segments;

  unsigned int allowed_missing = nlayers - rlayers;
  cout<<"allowed missing "<< allowed_missing<<endl;


  	// l is not actual layer number, it is ith layer
  	for (unsigned int l = 0; l < nlayers; ++l) {
    	layer_sorted[l].clear();
  	}
 
//  cout<<"track.hits.size "<<track.hits.size()<<endl;
  for (unsigned int i = 0; i < track.hits.size(); ++i) {
    Cluster3D hit = track.hits[i];
    unsigned int layer = (unsigned int) hit.get_layer();
    if (layer > (nlayers-1)) continue;
    if (!forward) layer = nlayers-layer-1;
    unsigned int min = (layer - allowed_missing);
    if (allowed_missing > layer) {
      min = 0;
    }
    for (unsigned int l = min; l <= layer; ++l) {
      layer_sorted[l].push_back(hit);
    }

  }

  for (unsigned int l = 0; l < nlayers; ++l) {
//    cout<<"layer_sorted["<<l<<"].size = "<< layer_sorted[l].size()<<endl;
    if (layer_sorted[l].size() == 0) {
      return 0;
    }
  }

  timeval t1, t2;
  double time1 = 0.;
  double time2 = 0.;

  gettimeofday(&t1, NULL);

  float ca_cos_ang_cut_diff = 1. - ca_cos_ang_cut;
  float ca_cos_ang_cut_diff_inv = 1. / ca_cos_ang_cut_diff;
  float ca_sin_ang_cut = sqrt(1. - ca_cos_ang_cut * ca_cos_ang_cut);

  std::vector<float> inv_layer;
  inv_layer.assign(nlayers, 1.);
  for (unsigned int l = 3; l < nlayers; ++l) {
    inv_layer[l] = 1. / (((float)l) - 2.);
  }

  float x1, x2, x3;
  float y1, y2, y3;
  float z1, z2, z3;
  float dx1, dx2, dx3;
  float dy1, dy2, dy3;
  float dz1, dz2, dz3;

  float kappa;
  float dkappa;

  float ux_mid;
  float uy_mid;
  float ux_end;
  float uy_end;

  float dzdl_1;
  float dzdl_2;
  float ddzdl_1;
  float ddzdl_2;


  float cur_kappa;
  float cur_dkappa;
  float cur_ux;
  float cur_uy;
  float cur_chi2;
  float chi2;

  unsigned int hit1;
  unsigned int hit2;
  unsigned int hit3;

  TrackSegment temp_segment;
  temp_segment.hits.assign(nlayers, 0);

  for (unsigned int i = 0; i < layer_sorted[0].size(); ++i) {
    for (unsigned int j = 0; j < layer_sorted[1].size(); ++j) {
      for (unsigned int k = 0; k < layer_sorted[2].size() ; ++k) {

	unsigned int layer0 = layer_sorted[0][i].get_layer(); 
	unsigned int layer1 = layer_sorted[1][j].get_layer();
	unsigned int layer2 = layer_sorted[2][k].get_layer();
	if (!forward) 
	{
	layer0 = nlayers-layer0-1;
	layer1 = nlayers-layer1-1;
	layer2 = nlayers-layer2-1;
	}
        if ((layer0 >= layer1) || (layer1 >= layer2)) {
          continue;
        }

        x1 = layer_sorted[0][i].get_x();
        y1 = layer_sorted[0][i].get_y();
        z1 = layer_sorted[0][i].get_z();

	// sigma ?= half pictch
        dx1 = 0.5*sqrt(12.0)*sqrt(layer_sorted[0][i].get_size(0,0));
        dy1 = 0.5*sqrt(12.0)*sqrt(layer_sorted[0][i].get_size(1,1));
        dz1 = 0.5*sqrt(12.0)*sqrt(layer_sorted[0][i].get_size(2,2));

        x2 = layer_sorted[1][j].get_x();
        y2 = layer_sorted[1][j].get_y();
        z2 = layer_sorted[1][j].get_z();

        dx2 = 0.5*sqrt(12.0)*sqrt(layer_sorted[1][j].get_size(0,0));
        dy2 = 0.5*sqrt(12.0)*sqrt(layer_sorted[1][j].get_size(1,1));
        dz2 = 0.5*sqrt(12.0)*sqrt(layer_sorted[1][j].get_size(2,2));

        x3 = layer_sorted[2][k].get_x();
        y3 = layer_sorted[2][k].get_y();
        z3 = layer_sorted[2][k].get_z();
        dx3 = 0.5*sqrt(12.0)*sqrt(layer_sorted[2][k].get_size(0,0));
        dy3 = 0.5*sqrt(12.0)*sqrt(layer_sorted[2][k].get_size(1,1));
        dz3 = 0.5*sqrt(12.0)*sqrt(layer_sorted[2][k].get_size(2,2));

	// layer of hit
        hit1 = i;
        hit2 = j;
        hit3 = k;

	calculate_kappa_tangents(x1, y1, z1, x2, y2, z2, x3, y3,
				z3, dx1, dy1, dz1, dx2, dy2, dz2,
				dx3, dy3, dz3, kappa, dkappa,
				ux_mid, uy_mid, ux_end, uy_end,
				dzdl_1, dzdl_2, ddzdl_1, ddzdl_2);

#ifdef _DEBUG_
	cout<<"Triplet : "<<endl;
	cout<<"kappa "<<kappa<< " dkappa "<<dkappa<<" ux_mid "<<ux_mid<<" uy_mid "<<" ux_end "<<ux_end<<
	" uy_end " <<uy_end<<" dzdl_1 "<<dzdl_1<<" dzdl_2 "<<dzdl_2<<" ddzdl_1 "<<ddzdl_1<<" ddzdl_2 "<<ddzdl_2<<endl;
#endif

          temp_segment.chi2 = pow(
              (dzdl_1 - dzdl_2) /
              (ddzdl_1 + ddzdl_2 + fabs(dzdl_1 * ca_sin_ang_cut)),2);
          if (temp_segment.chi2 > ca_chi2_layer_cut) continue;
          temp_segment.ux = ux_end;
          temp_segment.uy = uy_end;
          temp_segment.kappa = kappa;
//          if (temp_segment.kappa > _hough_space->get_kappa_max()) continue;
          temp_segment.dkappa = dkappa;
          temp_segment.hits[0] = hit1;
          temp_segment.hits[1] = hit2;
          temp_segment.hits[2] = hit3;
          temp_segment.n_hits = 3;
          unsigned int outer_layer =
          	layer_sorted[2][temp_segment.hits[2]].get_layer();
	  if (!forward) outer_layer = nlayers-outer_layer-1;
	  // make sure we have required number of layers with hits
          if ((outer_layer - 2) > allowed_missing) continue;
	  // finish up if required number of layers is reached,
          if ((nlayers - 3) <= allowed_missing) {
          	complete_segments.push_back(temp_segment);
          }
          if (next_seg->size() == next_seg_size) { // first new segment
          	next_seg->push_back(temp_segment);
          	next_seg_size += 1;
          } else { // next new segments
          	(*next_seg)[next_seg_size] = temp_segment;
          	next_seg_size += 1;
          }
	
      }
    }
  }

  swap(cur_seg, next_seg);
  swap(cur_seg_size, next_seg_size);

  cout<<"number of segments from first 3 layers : " << cur_seg_size<<endl;
  unsigned int which_seg;

  // add hits to segments layer-by-layer, cutting out bad segments
  for (unsigned int l = 3; l < nlayers; ++l) {
	if (l == (nlayers - 1)) {
//	ca_chi2_cut_layer = 0.25* ca_chi2_cut;// 2.*0.25 = 0.8 less loose cut after adding all clusters 
	}
    	next_seg_size = 0;
    for (unsigned int i = 0; i < cur_seg_size; ++i) {
      for (unsigned int j = 0; j < layer_sorted[l].size(); ++j) {

	unsigned int layer0 = layer_sorted[l - 1][(*cur_seg)[i].hits[l - 1]].get_layer();
	unsigned int layer1 = layer_sorted[l][j].get_layer();
	if (!forward) 
	{
		layer0 = nlayers-layer0-1;
		layer1 = nlayers-layer1-1;
	}
	if (layer0 >= layer1) continue;

        x1 = layer_sorted[l - 2][(*cur_seg)[i].hits[l - 2]].get_x();
       	y1 = layer_sorted[l - 2][(*cur_seg)[i].hits[l - 2]].get_y();
        z1 = layer_sorted[l - 2][(*cur_seg)[i].hits[l - 2]].get_z();
        x2 = layer_sorted[l - 1][(*cur_seg)[i].hits[l - 1]].get_x();
        y2 = layer_sorted[l - 1][(*cur_seg)[i].hits[l - 1]].get_y();
        z2 = layer_sorted[l - 1][(*cur_seg)[i].hits[l - 1]].get_z();
        x3 = layer_sorted[l][j].get_x();
        y3 = layer_sorted[l][j].get_y();
        z3 = layer_sorted[l][j].get_z();

        dx1 = 0.5*sqrt(12.0)*sqrt(layer_sorted[l - 2][(*cur_seg)[i].hits[l - 2]].get_size(0,0));
        dy1 = 0.5*sqrt(12.0)*sqrt(layer_sorted[l - 2][(*cur_seg)[i].hits[l - 2]].get_size(1,1));
        dz1 = 0.5*sqrt(12.0)*sqrt(layer_sorted[l - 2][(*cur_seg)[i].hits[l - 2]].get_size(2,2));
        dx2 = 0.5*sqrt(12.0)*sqrt(layer_sorted[l - 1][(*cur_seg)[i].hits[l - 1]].get_size(0,0));
        dy2 = 0.5*sqrt(12.0)*sqrt(layer_sorted[l - 1][(*cur_seg)[i].hits[l - 1]].get_size(1,1));
        dz2 = 0.5*sqrt(12.0)*sqrt(layer_sorted[l - 1][(*cur_seg)[i].hits[l - 1]].get_size(2,2));
        dx3 = 0.5*sqrt(12.0)*sqrt(layer_sorted[l][j].get_size(0,0));
        dy3 = 0.5*sqrt(12.0)*sqrt(layer_sorted[l][j].get_size(1,1));
        dz3 = 0.5*sqrt(12.0)*sqrt(layer_sorted[l][j].get_size(2,2));

        cur_kappa = (*cur_seg)[i].kappa;
        cur_dkappa = (*cur_seg)[i].dkappa;
        cur_ux = (*cur_seg)[i].ux;
        cur_uy = (*cur_seg)[i].uy;
        cur_chi2 = (*cur_seg)[i].chi2;

        which_seg = i;
        hit1 = j;

	calculate_kappa_tangents(
              x1, y1, z1, x2, y2, z2, x3, y3, z3, 
	      dx1, dy1, dz1, dx2, dy2, dz2, dx3, dy3, dz3, 
	      kappa, dkappa, ux_mid, uy_mid, ux_end, uy_end, 
	      dzdl_1, dzdl_2, ddzdl_1, ddzdl_2, 
	      ca_sin_ang_cut, ca_cos_ang_cut_diff_inv,
              cur_kappa, cur_dkappa, cur_ux, cur_uy, cur_chi2, chi2);

#ifdef _DEBUG_
	cout<<"Extended layers for segment "<<which_seg<<endl;
	cout<<"kappa "<<kappa<< " dkappa "<<dkappa
	<<" ux_mid "<<ux_mid<<" uy_mid "<<" ux_end "<<ux_end<<" uy_end " <<uy_end
	<<" dzdl_1 "<<dzdl_1<<" dzdl_2 "<<dzdl_2<<" ddzdl_1 "<<ddzdl_1<<" ddzdl_2 "<<ddzdl_2
	<<" cur_chi2 "<<cur_chi2<<" chi2 "<<chi2<<" chi2*inv_layer "<<l <<" "<<chi2*inv_layer[l]
	<<" chi2_cut" <<ca_chi2_layer_cut<<endl;
#endif
	if (chi2 * inv_layer[l] < ca_chi2_layer_cut) {
              temp_segment.chi2 = chi2;
              temp_segment.ux = ux_end;
              temp_segment.uy = uy_end;
              temp_segment.kappa = kappa;
              //if (temp_segment.kappa > _hough_space->get_kappa_max()) continue;
              temp_segment.dkappa = dkappa;
              for (unsigned int ll = 0; ll < l; ++ll) {
                temp_segment.hits[ll] = (*cur_seg)[which_seg].hits[ll];
              }
              temp_segment.hits[l] = hit1;
              unsigned int outer_layer =
                  layer_sorted[l][temp_segment.hits[l]].get_layer();
	      if (!forward) outer_layer = nlayers-outer_layer-1;
              temp_segment.n_hits = l + 1;
	      // finish up if required number of layers is reached
              if ((nlayers - (l + 1)) <= allowed_missing) {
                complete_segments.push_back(temp_segment);
              }
	      // make sure we have required number of layers with hits
              if ((outer_layer - l) > allowed_missing) {
                continue;
              }
              if (next_seg->size() == next_seg_size) { // first new segment
                next_seg->push_back(temp_segment);
                next_seg_size += 1;
              } else { // next new segments
                (*next_seg)[next_seg_size] = temp_segment;
                next_seg_size += 1;
              }
         }
      }// j 
    }// i

    swap(cur_seg, next_seg);
    swap(cur_seg_size, next_seg_size);

  }//l = 3

  cout<<"number of complete segments : " << complete_segments.size()<<endl;
  cout<<"number of current segments : "<< cur_seg_size<<endl;

  // copy complete segments over to current segments
  for (unsigned int i = 0; i < complete_segments.size(); ++i) {
    if (cur_seg->size() == cur_seg_size) {
      cur_seg->push_back(complete_segments[i]);
      ++cur_seg_size;
    } else {
      (*cur_seg)[cur_seg_size] = complete_segments[i];
      ++cur_seg_size;
    }
  }

  gettimeofday(&t2, NULL);
  time1 = ((double)(t1.tv_sec) + (double)(t1.tv_usec) / 1000000.);
  time2 = ((double)(t2.tv_sec) + (double)(t2.tv_usec) / 1000000.);
  CAtime += (time2 - time1);

  std::set<unsigned int> comp1;
  std::set<unsigned int> comp2;

  cout<<"number of segments generated "<< cur_seg_size<<endl;
  if (cur_seg_size==0  || cur_seg_size>10000) return 1;
  for (unsigned int i = cur_seg_size-1; i>0; --i){
    if ((*cur_seg)[i].n_hits==0) continue;
    comp1.clear();
    temp_combo.assign((*cur_seg)[i].n_hits, 0);
    for (unsigned int l = 0; l < (*cur_seg)[i].n_hits; ++l) {
    temp_combo[l] = layer_sorted[l][(*cur_seg)[i].hits[l]].get_id();
    comp1.insert(temp_combo[l]);
    }

    sort(temp_combo.begin(), temp_combo.end());
    set<vector<unsigned int> >::iterator it = combos.find(temp_combo);
    if (it != combos.end()) {
      (*cur_seg)[i].n_hits = 0.;
    }

//    if (combos.size() > 10000) {
    if (combos.size() > 100000) { 
      combos.clear();
    }
    combos.insert(temp_combo);

    for (unsigned int j = i-1; j>=0; --j){
      comp2.clear();
      for (unsigned int m = 0; m < (*cur_seg)[j].n_hits; ++m){
      comp2.insert(layer_sorted[m][(*cur_seg)[j].hits[m]].get_id());
      };
      for (std::set<unsigned int>::iterator it=comp1.begin(); it!=comp1.end(); ++it){
      auto it2 = comp2.find(*it);
      if (it2 != comp2.end()) comp2.erase(*it2);
      }	
      if (comp2.empty()) {
	(*cur_seg)[j].n_hits = 0;
      }
      if (j==0) break;
    }
  }

  unsigned int nsegs = cur_seg_size;
  for (unsigned int i = 0; i<cur_seg_size; i++) {
  if ((*cur_seg)[i].n_hits ==0) --nsegs;
  }
  cout<<"number of segments to be processed "<< nsegs<<endl;

  Track3D temp_track;
  temp_track.hits.assign(nlayers, Cluster3D());

  std::vector<Track3D> best_track;
  std::vector<HelixTrackState> best_track_state;
  float best_chi2 = 9999;
  for (unsigned int i = 0; i< cur_seg_size; ++i) {

    if ((*cur_seg)[i].n_hits ==0) continue;

#ifdef _DEBUG_
      cout<<"segment " <<i <<endl;
#endif
      temp_track.hits.assign((*cur_seg)[i].n_hits, Cluster3D());

      for (unsigned int l = 0; l < (*cur_seg)[i].n_hits; ++l) {
        temp_track.hits[l] = layer_sorted[l][(*cur_seg)[i].hits[l]];
      }

      unsigned int ninner_hits =0;
      for (unsigned int n = 0; n < temp_track.hits.size(); ++n){
        if (temp_track.hits[n].get_layer()<3)
	++ninner_hits;
      }
     
      if(require_inner_hits && ninner_hits<2) continue;

      gettimeofday(&t1, NULL);

      float init_chi2 = temp_track.fit_track();
#ifdef _DEBUG_
      cout<<"chi2 from fit_track "<<init_chi2
  	  <<" kappa " <<temp_track.kappa <<" phi "<<temp_track.phi<<" d "<<temp_track.d
	  <<" z0 "<<temp_track.z0<<" dzdl "<<temp_track.dzdl<< endl;
#endif

    // not being used for the time being
    if (init_chi2 > fast_chi2_cut_max) { 
      if (init_chi2 > fast_chi2_cut_par0 +
                          fast_chi2_cut_par1 / kappa_to_pt(temp_track.kappa)) {
        gettimeofday(&t2, NULL);
        time1 = ((double)(t1.tv_sec) + (double)(t1.tv_usec) / 1000000.);
        time2 = ((double)(t2.tv_sec) + (double)(t2.tv_usec) / 1000000.);
        KALtime += (time2 - time1);
        continue;
      }
    }

    HelixTrackState state;
    state.phi = temp_track.phi;
    if (state.phi < 0.) {
      state.phi += 2. * M_PI;
    }
    state.d = temp_track.d;
    state.kappa = temp_track.kappa;
    state.nu = sqrt(state.kappa);
    state.z0 = temp_track.z0;
    state.dzdl = temp_track.dzdl;
    state.C = Matrix<float, 5, 5>::Zero(5, 5);
    state.C(0, 0) = pow(0.01, 2.);
    state.C(1, 1) = pow(0.01, 2.);
    state.C(2, 2) = pow(0.01 * state.nu, 2.);
    state.C(3, 3) = pow(0.05, 2.);
    state.C(4, 4) = pow(0.05, 2.);
    state.chi2 = 0.;
    state.position = 0;
    state.x_int = 0.;
    state.y_int = 0.;
    state.z_int = 0.;
 
    unsigned int nfits = 0;
    for (unsigned int h = 0; h < temp_track.hits.size(); ++h) {
      _kalman->addHit(temp_track.hits[h], state);
      nfits += 1;
    }
    cout<<"z0 after kalman "<<state.z0<<endl;

    // fudge factor for non-gaussian hit sizes
    state.C *= 3.;
    state.chi2 *= 6.;

    gettimeofday(&t2, NULL);
    time1 = ((double)(t1.tv_sec) + (double)(t1.tv_usec) / 1000000.);
    time2 = ((double)(t2.tv_sec) + (double)(t2.tv_usec) / 1000000.);
    KALtime += (time2 - time1);

    if (!(temp_track.kappa == temp_track.kappa)) {
      continue;
    }
/*
    if (temp_track.kappa > _hough_space->get_kappa_max()) {
      continue;
    }
*/
    if (!(state.chi2 == state.chi2)) {
      continue;
    }
    if (state.chi2 / (2. * ((float)(temp_track.hits.size())) - 5.) > ca_chi2_cut) {
      continue;
    }
/*
    if (cut_on_dca == true) {
      if (fabs(temp_track.d) > dca_cut) {
        continue;
      }
      if (fabs(temp_track.z0) > dca_cut) {
        continue;
      }
    }
*/
    if (best_chi2 > state.chi2){
      if (!best_track.empty()){
      best_track.pop_back();
      best_track_state.pop_back();
      }
      best_track.push_back(temp_track);
      best_track_state.push_back(state);
    }
  }
  
  if (best_track.empty()) return 1;

  ca_tracks.push_back(best_track.back());
  ca_track_states.push_back(best_track_state.back());
  cout <<"ca track added, chi2 =  "<< (best_track_state.back().chi2)/(2. * ((float)(temp_track.hits.size())) - 5.) <<" z0 = "<<best_track_state.back().z0<<endl;


//    if ((remove_hits == true) && (state.chi2 < chi2_removal_cut) &&
//        (temp_track.hits.size() >= n_removal_hits)) {
	temp_track = best_track.back();

	if  (remove_hits){
	for (unsigned int i = 0; i < temp_track.hits.size(); ++i) {
		if (!remove_inner_hits && temp_track.hits[i].get_layer()<3) continue;
		auto search = _hits_used.find(temp_track.hits[i].get_id());
		if(search != _hits_used.end())
        	{
        	_hits_used.find(temp_track.hits[i].get_id())->second = true;
        	}
	}
	}

        segments1.clear();
        segments2.clear();


  return 1;
}

int CellularAutomaton_v1::extend_single_triplet(Track3D& track){
	/***
	 if INTT is used, compute the projected z position of an INTT cluster utilizing the kappa 
	 of current segment consisting of triplets (z0,z1,z2), and then update the kappa of new quadraplet.  
          Input: triplets in Track3D, Output: z_proj (projected z3) 

         all of the helical parameters have to be re-calculated ====> Save Track3D on node!!         
         store z_proj in proj_hits_map as z (z3=z_proj) and in hit3d as z in current segment. 
	***/
	if (verbosity>0) cout<<"entered extend single track"<<endl;
	// sort out clusters by layer, triplets in the first 3 layers
	for(int i=0; i<3; ++i) layer_sorted[i].clear();
	layer_sorted[0].push_back(track.hits[0]);
	layer_sorted[1].push_back(track.hits[1]);
	layer_sorted[2].push_back(track.hits[2]);
	int trackid = track.get_id();
/*
        for (std::map<unsigned int,Cluster3D>::iterator jt = _hits_map.begin();
                	jt!= _hits_map.end();
                        ++jt) {

                Cluster3D hit = jt->second;
                layer_sorted[l].push_back(hit);
                }
 */
        for (unsigned int l = 0; l< (nlayers+3); ++l){
        if(verbosity>1)cout<<"layer_sorted["<<l<<"].size = "<< layer_sorted[l].size()<<endl;
                if (layer_sorted[l].size() == 0) {
                return 0;
                }
        }

        // a triplet stored in a segment can be used multiple times as layers are added
        unsigned int cur_seg_size = 0;
        unsigned int next_seg_size = 0;
	std::map<unsigned int, TrackSegment> cur_seg;
	std::map<unsigned int, TrackSegment> next_seg;
        std::map<unsigned int, Cluster3D> _proj_hits_map;
	std::map<unsigned int, unsigned int> _proj_orig_map;
	// this will be needed when allowing missing hits
//        std::vector<TrackSegment> complete_segments;

        float ca_cos_ang_cut_diff = 1. - ca_cos_ang_cut;
        float ca_cos_ang_cut_diff_inv = 1. / ca_cos_ang_cut_diff;
        float ca_sin_ang_cut = sqrt(1. - ca_cos_ang_cut * ca_cos_ang_cut);

        std::vector<float> inv_layer;
        inv_layer.assign(nlayers, 1.);
//        for (unsigned int l = 3; l < nlayers+3; ++l) {
//        inv_layer[l] = 1. / (((float)(l)) - 2.);
//        }
	inv_layer[3]=1.;

        float x1, x2, x3;
        float y1, y2, y3;
        float z1, z2, z3;
        float dx1, dx2, dx3;
        float dy1, dy2, dy3;
        float dz1, dz2, dz3;

        float kappa;
        float dkappa;

        float ux_mid;
        float uy_mid;
        float ux_end;
        float uy_end;

        float dzdl_1;
        float dzdl_2;
        float ddzdl_1;
        float ddzdl_2;

#ifdef _FULL_TEST_
        float cur_kappa;
        float cur_dkappa;
        float cur_ux;
        float cur_uy;
        float cur_chi2;
        float chi2;
#endif

        unsigned int hit1;
        unsigned int hit2;
        unsigned int hit3;

        TrackSegment temp_segment;
        temp_segment.hits.assign(nlayers+3, 0);
//	temp_segment.projzs.assign(nlayers+3,-999.);
	std::vector<float> projzs;
	projzs.assign(nlayers+3,-999.);

	x1 = track.hits[0].get_x();
	y1 = track.hits[0].get_y();
	z1 = track.hits[0].get_z();
	dx1 = track.hits[0].get_size(0,0);
	dy1 = track.hits[0].get_size(1,1);
	dz1 = track.hits[0].get_size(2,2);

	x2 = track.hits[1].get_x();
	y2 = track.hits[1].get_y();
	z2 = track.hits[1].get_z();
	dx2 = track.hits[1].get_size(0,0);
	dy2 = track.hits[1].get_size(1,1);
	dz2 = track.hits[1].get_size(2,2);

	x3 = track.hits[2].get_x();
	y3 = track.hits[2].get_y();
	z3 = track.hits[2].get_z();
	dx3 = track.hits[2].get_size(0,0);
	dy3 = track.hits[2].get_size(1,1);
	dz3 = track.hits[2].get_size(2,2);

	hit1 = 0; hit2=0; hit3=0; // hit index in layer_sorted[layer][hit_index], saved in segments

        calculate_kappa_tangents(x1, y1, z1, x2, y2, z2, x3, y3,
                                z3, dx1, dy1, dz1, dx2, dy2, dz2,
                                dx3, dy3, dz3, kappa, dkappa,
                                ux_mid, uy_mid, ux_end, uy_end,
                                dzdl_1, dzdl_2, ddzdl_1, ddzdl_2);

#ifdef _DEBUG_
        cout<<"Triplet : "<<endl;
        cout<<"kappa "<<kappa<< " dkappa "<<dkappa<<" ux_mid "<<ux_mid<<" uy_mid "<<" ux_end "<<ux_end<<
        " uy_end " <<uy_end<<" dzdl_1 "<<dzdl_1<<" dzdl_2 "<<dzdl_2<<" ddzdl_1 "<<ddzdl_1<<" ddzdl_2 "<<ddzdl_2<<endl;
#endif

        temp_segment.chi2 = pow(
        	(dzdl_1 - dzdl_2) /
                (ddzdl_1 + ddzdl_2 + fabs(dzdl_1 * ca_sin_ang_cut)),2);
        temp_segment.ux = ux_end;
        temp_segment.uy = uy_end;
        temp_segment.kappa = kappa;
        temp_segment.dkappa = dkappa;

        temp_segment.hits[0] = hit1;// hit index in layer_sorted for each layer.  
        temp_segment.hits[1] = hit2;
        temp_segment.hits[2] = hit3;
        temp_segment.n_hits = 3;

//        if (cur_seg.size() == cur_seg_size) {
        cur_seg.insert(make_pair(cur_seg_size,temp_segment));
        ++cur_seg_size;
//        } else {
//        (cur_seg)[cur_seg_size] = temp_segment;
//        ++cur_seg_size;
//        }

	//set up Kalman Filter class and get helix parameters for a triplet to be used to get z projection
        float init_chi2 = track.fit_track();
        if (verbosity>3) cout<<"chi2 from fit_track "<<init_chi2 <<" kappa "<< track.kappa <<endl;
		// get helix parameters to be used for kalman filter 
                HelixTrackState state;
                state.phi = track.phi;
                if (state.phi < 0.) {
                state.phi += 2. * M_PI;
                }
                state.d = track.d;
                state.kappa = track.kappa;
                state.nu = sqrt(state.kappa);
                state.z0 = track.z0;
                state.dzdl = track.dzdl;
                state.C = Matrix<float, 5, 5>::Zero(5, 5);
                state.C(0, 0) = pow(0.01, 2.);
                state.C(1, 1) = pow(0.01, 2.);
                state.C(2, 2) = pow(0.01 * state.nu, 2.);
                state.C(3, 3) = pow(0.05, 2.);
                state.C(4, 4) = pow(0.05, 2.);
                state.chi2 = 0.;
                state.position = 0;
                state.x_int = 0.;
                state.y_int = 0.;
                state.z_int = 0.;

                unsigned int nfits = 0;
                for (unsigned int h = 0; h < track.hits.size(); ++h) {
                _kalman->addHit(track.hits[h], state);
                nfits += 1;
		}
//                state.C *= 3.;
//                state.chi2 *= 6.;

	// now we have initial helical parameters saved in state 
	// go to next layer and get z projection   
        unsigned int added_next_segments = 0;
	float best_chi2=999.;
	Track3D temp_track;
	Track3D best_track;
	for(unsigned int l = 3; l<(nlayers+3); ++l ){
		if (verbosity>1)cout<<"extend_triplet : layer "<<l<<" type "<<_ladder_type[(l-3)*2]<<" (z seg : -1, phi seg : -2)" <<endl;
                next_seg_size = 0;
                next_seg.clear();
                // Loop over current segments
                        for (unsigned int i = 0; i < cur_seg_size; ++i) {
                        if ((cur_seg)[i].n_hits ==0) continue;
                        added_next_segments = 0;

		// first 2 clusters from previous 
                        // get positions  from layer_sorted 
                        x1 = layer_sorted[l - 2][(cur_seg)[i].hits[l - 2]].get_x();
                        y1 = layer_sorted[l - 2][(cur_seg)[i].hits[l - 2]].get_y();
                        z1 = layer_sorted[l - 2][(cur_seg)[i].hits[l - 2]].get_z();
			if ( use_projection && ((l-2)>=3)){
 			   if(_ladder_type[(l-5)*2]==-2){
				z1= _proj_hits_map[(cur_seg)[i].hits[l-2]].get_z();
			        //z1 = (cur_seg)[i].proj[l-5];
			   } else {
				x1 = _proj_hits_map[(cur_seg[i].hits[l-2])].get_x();
				x2 = _proj_hits_map[(cur_seg[i].hits[l-2])].get_y(); 
				//float phi_proj = (cur_seg)[i].proj[l-5];
			        //float r = sqrt(x1*x1+y1*y1);
                                //x1 = r*cos(phi_proj); y1 = r*sin(phi_proj);
			   }
			}			
                        dx1 = 0.5*sqrt(12.0)*sqrt(layer_sorted[l - 2][(cur_seg)[i].hits[l - 2]].get_size(0,0));
                        dy1 = 0.5*sqrt(12.0)*sqrt(layer_sorted[l - 2][(cur_seg)[i].hits[l - 2]].get_size(1,1));
                        dz1 = 0.5*sqrt(12.0)*sqrt(layer_sorted[l - 2][(cur_seg)[i].hits[l - 2]].get_size(2,2));
                        // get positions from hits_map
//                        auto searchl2 = _hits_map.find((cur_seg)[i].hits[l-2]);
//                        Cluster3D clusterl2 = searchl2->second;
//                        x1 = clusterl2.get_x();
//                        y1 = clusterl2.get_y();
//                        z1 = clusterl2.get_z();
//                        dx1 = 0.5*sqrt(12.0)*sqrt(clusterl2.get_size(0,0));
//                        dy1 = 0.5*sqrt(12.0)*sqrt(clusterl2.get_size(1,1));
//                        dz1 = 0.5*sqrt(12.0)*sqrt(clusterl2.get_size(2,2));
                        x2 = layer_sorted[l - 1][(cur_seg)[i].hits[l - 1]].get_x();
                        y2 = layer_sorted[l - 1][(cur_seg)[i].hits[l - 1]].get_y();
                        z2 = layer_sorted[l - 1][(cur_seg)[i].hits[l - 1]].get_z();
                        if ( use_projection && ((l-1)>=3)){
                           if(_ladder_type[(l-4)*2]==-2){
				z2= _proj_hits_map[(cur_seg)[i].hits[l-1]].get_z();
                           } else {
                                x2 = _proj_hits_map[(cur_seg[i].hits[l-1])].get_x();
                                y2 = _proj_hits_map[(cur_seg[i].hits[l-1])].get_y();
                                //float phi_proj = (cur_seg)[i].proj[l-4];
                                //float r = sqrt(x2*x2+y2*y2);
                                //x2 = r*cos(phi_proj); y2 = r*sin(phi_proj);
                           }
                        }
                        dx2 = 0.5*sqrt(12.0)*sqrt(layer_sorted[l - 1][(cur_seg)[i].hits[l - 1]].get_size(0,0));
                        dy2 = 0.5*sqrt(12.0)*sqrt(layer_sorted[l - 1][(cur_seg)[i].hits[l - 1]].get_size(1,1));
                        dz2 = 0.5*sqrt(12.0)*sqrt(layer_sorted[l - 1][(cur_seg)[i].hits[l - 1]].get_size(2,2));

                	float z_proj = -999.;//_kalman->get_z_projection(l);
			float phi_proj = -999;
			best_chi2 = 999.;
				// loop over clusters on current layer
				for(unsigned int ct = 0; ct < layer_sorted[l].size(); ++ct){
                                Cluster3D hit3d = layer_sorted[l][ct];
                                hit1 = ct; // hit index in layer_sort

                                x3 = hit3d.get_x();
                                y3 = hit3d.get_y();
				z3 = hit3d.get_z();
                                float phi_cur = shift_phi_range(atan2(y3,x3));

				// if a cluster falls within the search window
				// save proj_z in the segment rather than Cluster3D
				z_proj = _kalman->getZProjection(hit3d, state);
				phi_proj= _kalman->getPhiProjection(hit3d,state);

				phi_proj = shift_phi_range(phi_proj);
				float phi_diff = phi_cur - phi_proj;
				if (phi_cur< M_PI/2. && phi_proj > 3*M_PI/2.) phi_diff += 2.*M_PI;
				else if (phi_cur >3*M_PI/2. && phi_proj<M_PI/2.) phi_diff -= 2.*M_PI;
				
				float zcut=0.;float phicut=0.;  
                                if (_ladder_type[(l-3)*2]==-2) { //phi segmentation
					zcut = ca_zproj_z_cut; phicut = ca_zproj_phi_cut; 
					if(l==4) phicut *=4.2;// multiple scattering
				} else if (_ladder_type[(l-3)*2]==-1){// z segmentation
					zcut = ca_pproj_z_cut; phicut = ca_pproj_phi_cut;
				}
					if(l==3){ var=phi_diff;}// z_proj-z3; 
	 	
				if(verbosity>1) cout<<"  z3 proj "<<z_proj<< " z3 "<< z3<< " phi proj "<<phi_proj<<" phi "<<phi_cur <<endl;
				if(verbosity>3) cout<<"    phi cut "<<phicut<<" phi diff "<<phi_diff<<" zcut "<<zcut<< " z diff "<<fabs(z3-z_proj)<<endl;
				if((fabs(phi_diff)> phicut) || (fabs(z3-z_proj)> zcut )) continue;
				
				if (use_projection && _ladder_type[(l-3)*2]== -2 )//PHG4SiliconTrackerDefs::SEGMENTATION_PHI = -2 
				{
					z3 = z_proj;
					hit2 = _proj_hits_map.size();
					Cluster3D projhit3d = hit3d;
                                        projhit3d.set_id(hit2);
                                        projhit3d.set_x(x3);
                                        projhit3d.set_y(y3);
                                        projhit3d.set_z(z3);
					_proj_hits_map.insert(make_pair(hit2, projhit3d));
					_proj_orig_map.insert(make_pair(hit2, hit1));
					hit1 = hit2;
					if (verbosity>1) cout<<"    using projection z"<<endl;
				}
                                if (use_projection && _ladder_type[(l-3)*2]==-1)//PHG4SiliconTrackerDefs::SEGMENTATION_Z = -1
                                {
					if(verbosity>1) cout<<"    using projection phi"<<endl;
					float r = sqrt(x3*x3+y3*y3);
					x3 = r*cos(phi_proj); y3 = r*sin(phi_proj);
                                        hit2 = _proj_hits_map.size();
                                        Cluster3D projhit3d =  hit3d;                
					projhit3d.set_id(hit2);
                                	projhit3d.set_x(x3);
                			projhit3d.set_y(y3);
                			projhit3d.set_z(z3);
                                        _proj_hits_map.insert(make_pair(hit2, projhit3d));
					_proj_orig_map.insert(make_pair(hit2, hit1));
					hit1 = hit2;
				}	
				dx3 = 0.5*sqrt(12.0)*sqrt(hit3d.get_size(0,0));
                                dy3 = 0.5*sqrt(12.0)*sqrt(hit3d.get_size(1,1));
                                dz3 = 0.5*sqrt(12.0)*sqrt(hit3d.get_size(2,2));
				//if (verbosity>2) cout<<"dz3 "<<dz3<< " distance "<<endl;
                                cur_kappa = (cur_seg)[i].kappa;
                                cur_dkappa = (cur_seg)[i].dkappa;
                                cur_ux = (cur_seg)[i].ux;
                                cur_uy = (cur_seg)[i].uy;
                                cur_chi2 = (cur_seg)[i].chi2;

				// update kappa, dzdl and chi2
                                chi2=9999;
                                calculate_kappa_tangents(
                                        x1, y1, z1, x2, y2, z2, x3, y3, z3,
                                        dx1, dy1, dz1, dx2, dy2, dz2, dx3, dy3, dz3,
                                        kappa, dkappa, ux_mid, uy_mid, ux_end, uy_end,
                                        dzdl_1, dzdl_2, ddzdl_1, ddzdl_2,
                                        ca_sin_ang_cut, ca_cos_ang_cut_diff_inv,
                                        cur_kappa, cur_dkappa, cur_ux, cur_uy, cur_chi2, chi2);

				if (verbosity>1) cout<<"     chi2 of quadraplet candidate for this triplet "<<chi2<<endl;
				// if passing chi2 from kappa_tangents, fit track
                                if (chi2 * inv_layer[l] < ca_chi2_layer_cut ) {
					// get chi2 from track fitting 
					temp_track.hits.assign(l+1,Cluster3D());
                                        for (unsigned int ll = 0; ll < l; ++ll) {
					  if(ll>=3 && use_projection){
					  temp_track.hits[ll] = layer_sorted[ll][_proj_orig_map[cur_seg[i].hits[ll]]];
					  }else{
                                          temp_track.hits[ll] = layer_sorted[ll][(cur_seg)[i].hits[ll]];
					  }
                                        }
					temp_track.hits[l]=hit3d;
					float temp_chi2 = temp_track.fit_track();
					if (verbosity>1) cout<<"     chi2 of temp_segment "<<temp_chi2<<endl;
					// add best hit to temp_segment
					if(best_chi2> chi2){
					best_chi2 = chi2;
					best_track = temp_track;
                                        temp_segment.chi2 = chi2;
                                        temp_segment.ux = ux_end;
                                        temp_segment.uy = uy_end;
                                        temp_segment.kappa = kappa;
                                        temp_segment.dkappa = dkappa;
                                        for (unsigned int ll = 0; ll < l; ++ll) {
                                        temp_segment.hits[ll] = (cur_seg)[i].hits[ll];
                                        }
                                        temp_segment.hits[l] = hit1;
                                        temp_segment.n_hits = l + 1;
					//temp_segment.proj[l-3]=proj;
					}

				} // chi2 cut
                                }// clusters on layer l
				if (verbosity>1) cout<<"best chi2 for this triplet"<< best_chi2<<endl;
					// add temp_segment(with best hits) to next segment
                                        if (next_seg.size() == next_seg_size) { // first new segment
                                                next_seg.insert(make_pair(next_seg_size,temp_segment));
                                                next_seg_size += 1;
                                        } else { // next new segments
                                        #ifdef _DEBUG_
                                        cout<<"Next segment size inconsistent "<<endl;
                                        #endif
                                                (next_seg)[next_seg_size] = temp_segment;
                                                next_seg_size += 1;
                                        }
                                        ++added_next_segments;
#ifdef _DEBUG_                                  
                                        cout<<"segment "<< i<< " added segment "<<added_next_segments<<endl
#endif

			}// current segments
		// select best hit and add to kalman filter to get updated z-projection
		// alternative: accumulate segments and until end of layers and add hits at a time
		// 		this results in repetative computation with kalman filter.   	
		if(!use_projection){
                _kalman->addHit(layer_sorted[l][temp_segment.hits[l]], state);
		}else{
		_kalman->addHit(_proj_hits_map[temp_segment.hits[l]], state);
		}
                cur_seg.swap(next_seg);
                swap(cur_seg_size, next_seg_size);
	}//nlayers


	// update track and state and then save
//	state.C *= 3.;
//	state.chi2 *= 6.;
 	//**** get chi2 of last layer******//
	//state.chi2 = best_chi2;// remove when tracking
  	//*********************************//

        // when extending to 2nd intt layer, create a list of best_track associated with each segment and get the best                
	best_track.set_id(trackid);
	ca_tracks.push_back(best_track);
        ca_track_states.push_back(state);
	if (verbosity>0) cout<<"extended ca track added, track id "<<best_track.get_id()<<endl;

/*
        std::vector<Track3D> best_track;
        std::vector<HelixTrackState> best_track_state;

        ca_tracks.push_back(best_track.back());
        ca_track_states.push_back(best_track_state.back());
*/
        cur_seg.clear();
        next_seg.clear();

  return 1;
}

int CellularAutomaton_v1::calculate_kappa_tangents(
			float x1, float y1, float z1, float x2, float y2, float z2, 
			float x3, float y3, float z3, 
			float dx1, float dy1, float dz1, float dx2, float dy2, float dz2,
                        float dx3, float dy3, float dz3, 
			float& kappa, float& dkappa,
                        float& ux_mid, float& uy_mid, float& ux_end, float& uy_end,
                        float& dzdl_1, float& dzdl_2, float& ddzdl_1, float& ddzdl_2)
{

	float D12 = sqrt(pow(x2-x1,2)+pow(y2-y1,2));
	float D23 = sqrt(pow(x3-x2,2)+pow(y3-y2,2));
	float D13 = sqrt(pow(x3-x1,2)+pow(y3-y1,2));
	kappa = 1./(D12*D23*D13);
	float num = (D12+D23+D13)*(D23+D13-D12)*(D12+D13-D23)*(D12+D23-D13);
	if (num<0) num = 0;
	num = sqrt(num);	
	kappa *=num;
	
	float kappa_inv = 1/kappa;
	float D12_inv = 1./D12;
	float D23_inv = 1./D23;
	float D13_inv = 1./D13;

	float dr1 = sqrt(pow(dx1,2)+pow(dy1,2));
	float dr2 = sqrt(pow(dx2,2)+pow(dy2,2));
	float dr3 = sqrt(pow(dx3,2)+pow(dy3,2));

	float dk1 = (dr1+dr2)* D12_inv*D12_inv;
	float dk2 = (dr2+dr3)* D23_inv*D23_inv;
	dkappa = dk1+dk2;

	float ux12 = (x2-x1)*D12_inv;
	float uy12 = (y2-y1)*D12_inv;
	float ux23 = (x3-x2)*D23_inv;
	float uy23 = (y3-y2)*D23_inv;
	float ux13 = (x3-x1)*D13_inv;
	float uy13 = (y3-y1)*D13_inv;

	// cos(alpha) = cos(alpha12 - alpha13)
	// sin(alpha) = sin(alpha12 - alpha13)
	float cosalpha = ux12*ux13 + uy12*uy13;
	float sinalpha = uy12*ux13 - ux12*uy13;
	// alpha23 + alpha  
	ux_mid = ux23 * cosalpha - uy23 * sinalpha;
	uy_mid = ux23 * sinalpha + uy23 * cosalpha;
	// alpha23 - alpha
	ux_end = ux23 * cosalpha + uy23 * sinalpha;
	uy_end = uy23 * cosalpha - ux23 * sinalpha;

	// dzdl = dz/sqrt(ds^2 + dz^2)
	float ds23 = 2.*kappa_inv*atan(sinalpha/(1.+ sqrt(1.-pow(sinalpha,2))));
	if (kappa<=0)  ds23 = D23;

	float dz23 = z3 - z2;
	dzdl_2 =  dz23/sqrt(pow(ds23,2) + pow(dz23,2));
	ddzdl_2 = (dz2 + dz3)*D23_inv;

	// sin(alpha) = sin(alpha13 -alpha23)
	sinalpha = ux13 *uy23 - ux23 *  uy13;
	float ds12 = 2.*kappa_inv*atan(sinalpha/(1.+sqrt(1.-pow(sinalpha,2))));
	if (kappa<=0) ds12 = D12;

	float dz12 = z2 - z1;
	dzdl_1 = dz12/sqrt(pow(ds12,2) + pow(dz12,2));
	ddzdl_1 = (dz1 + dz2) * D12_inv;

	return 1;

}

int CellularAutomaton_v1::calculate_kappa_tangents(
                        float x1, float y1, float z1, float x2, float y2, float z2,
                        float x3, float y3, float z3,
                        float dx1, float dy1, float dz1, float dx2, float dy2, float dz2,
                        float dx3, float dy3, float dz3,
                        float& kappa, float& dkappa,
                        float& ux_mid, float& uy_mid, float& ux_end, float& uy_end,
                        float& dzdl_1, float& dzdl_2, float& ddzdl_1, float& ddzdl_2,
			float ca_sin_ang_cut, float ca_cos_ang_cut_diff_inv,
              		float cur_kappa, float cur_dkappa, float cur_ux, float cur_uy, 
			float cur_chi2, float& chi2)
{


        float D12 = sqrt(pow(x2-x1,2)+pow(y2-y1,2));
        float D23 = sqrt(pow(x3-x2,2)+pow(y3-y2,2));
        float D13 = sqrt(pow(x3-x1,2)+pow(y3-y1,2));
        kappa = 1./(D12*D23*D13);
        float num = (D12+D23+D13)*(D23+D13-D12)*(D12+D13-D23)*(D12+D23-D13);
        if (num<0) num = 0;
        num = sqrt(num);
        kappa *=num;

        float kappa_inv = 1/kappa;
        float D12_inv = 1./D12;
        float D23_inv = 1./D23;
        float D13_inv = 1./D13;

        float dr1 = sqrt(pow(dx1,2)+pow(dy1,2));
        float dr2 = sqrt(pow(dx2,2)+pow(dy2,2));
        float dr3 = sqrt(pow(dx3,2)+pow(dy3,2));

        float dk1 = (dr1+dr2)* D12_inv*D12_inv;
        float dk2 = (dr2+dr3)* D23_inv*D23_inv;
        dkappa = dk1+dk2;

        float ux12 = (x2-x1)*D12_inv;
        float uy12 = (y2-y1)*D12_inv;
        float ux23 = (x3-x2)*D23_inv;
        float uy23 = (y3-y2)*D23_inv;
        float ux13 = (x3-x1)*D13_inv;
        float uy13 = (y3-y1)*D13_inv;

        // cos(alpha) = cos(alpha12 - alpha13)
        // sin(alpha) = sin(alpha12 - alpha13)
        float cosalpha = ux12*ux13 + uy12*uy13;
        float sinalpha = uy12*ux13 - ux12*uy13;
        // alpha23 + alpha  
        ux_mid = ux23 * cosalpha - uy23 * sinalpha;
        uy_mid = ux23 * sinalpha + uy23 * cosalpha;
        // alpha23 - alpha
        ux_end = ux23 * cosalpha + uy23 * sinalpha;
        uy_end = uy23 * cosalpha - ux23 * sinalpha;

        // dzdl = dz/sqrt(ds^2 + dz^2)
        float ds23 = 2.*kappa_inv*atan(sinalpha/(1.+ sqrt(1.-pow(sinalpha,2))));
        if (kappa<=0)  ds23 = D23;

        float dz23 = z3 - z2;
        dzdl_2 =  dz23/sqrt(pow(ds23,2) + pow(dz23,2));
        ddzdl_2 = (dz2 + dz3)*D23_inv;

        // sin(alpha) = sin(alpha13 -alpha23)
        sinalpha = ux13 *uy23 - ux23 *  uy13;
        float ds12 = 2.*kappa_inv*atan(sinalpha/(1.+sqrt(1.-pow(sinalpha,2))));
        if (kappa<=0) ds12 = D12;

        float dz12 = z2 - z1;
        dzdl_1 = dz12/sqrt(pow(ds12,2) + pow(dz12,2));
        ddzdl_1 = (dz1 + dz2) * D12_inv;

	float kappa_diff = cur_kappa - kappa;
	float n_dk = cur_dkappa	+ dkappa + ca_sin_ang_cut * kappa;
	float chi2_kappa = pow(kappa_diff,2)/pow(n_dk,2);

	float cos_scatter = cur_ux * ux_mid + cur_uy * uy_mid;
	float chi2_ang = pow((1-cos_scatter)*ca_cos_ang_cut_diff_inv,2);

	float sin_scatter = dzdl_1 * ca_sin_ang_cut;
	float chi2_dzdl = 0.5*pow((dzdl_1-dzdl_2)/(ddzdl_1+ddzdl_2+fabs(sin_scatter)),2);

	chi2 = cur_chi2 + chi2_ang + chi2_kappa + chi2_dzdl;
	return 1;
}


float CellularAutomaton_v1::kappa_to_pt(float kappa) {
        return _pt_rescale * _mag_field / 333.6 / kappa;
}

float CellularAutomaton_v1::shift_phi_range(float _phi){

        if (_phi < 0.) _phi += 2.*M_PI;
        return _phi;
}


