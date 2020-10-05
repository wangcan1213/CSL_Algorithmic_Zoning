/***
* Name: workshop
* Author: Can Wang
* Description: 
* Tags: Tag1, Tag2, TagN
***/

model workshop

/* Insert your model definition here */

global {
	file<geometry> grid_shp <- file<geometry>("../includes/grid80.shp");
	file<geometry> block_shp <- file<geometry>("../includes/greater_area.shp");
	file workers_csv <- csv_file("../includes/workers.csv", true);
	string blockgroup_distance_matrix_json_path <- "../includes/blockgroup_distance.json";
	string blockgroup_to_grid_distance_matrix_json_path <- "../includes/blockgroup_to_grid_distance.json";
	string grid_distance_matrix_json_path <- "../includes/grid_distance.json";
	geometry shape <- envelope(block_shp);
	geometry kendall_envelope;
	list<string> kendall_geoid_list <- ['250173526001', '250173521022', '250173521021', 
										'250173523002', '250173523001', '250173524001', 
										'250173524002', '250173526002', '250173522001'];
	map<string,float> population_discount <- ['250173522001'::0.5, '250173521022'::0.5, '250173521021'::0.2];
	map<string, map<string, float>> blockgroup_distance_matrix_map;
	map<string, map<string, float>> blockgroup_to_grid_distance_matrix_map;
	map<string, map<string, float>> grid_distance_matrix_map;
	map<string, list<string>> grid_neighbors_map;
	map<string,blockgroup> blockgroup_lookup_map;
	list<blockgroup> kendall_blockgroups;
	list<blockgroup> nonkendall_blockgroups;
	list<string> nonkendall_geoid_list;
	blockgroup kendall_virtual_block;
	float far_discount <- 0.5;    //it seems that current landuse plan provides too many residences
	list<map> move_stat <- []; 
	bool consider_nonkendall_resident_capacity <- false;
	int max_cycle <- 12;
	developer the_developer;
	float start_rent_coeff <- 0.95;
	
	int lower_nb_represented_persons <- 3;     //mainly used for kendall landuse, and as the unit for a real move
	int upper_nb_represented_persons <- 120;   //mainly used for outside blockgroup
	int nb_nonkendall_blocks_in_hom_loc_choice <- 8;
	float ratio_of_people_considering_move_at_each_step <- 0.05;
	float mean_monthly_rent_low_income <- 1112.0; //25% percentile for all rents, used as rent of current residence for low inc person if block mean rent unavailabe
	float mean_monthly_rent_high_income <- 1429.0; //75% percentile for all rents, used as rent of current residence for high inc person if block mean rent unavailabe
	float beta_work_allocation <- -0.5;
	float beta_work_allocation_for_kendall_residents <- -0.6;
	float residence_area_small_scale <- 40.0;
	float residence_area_large_scale <- 80.0;
	float office_area_small_scale <- 20.0;
	float office_area_large_scale <- 40.0;
	float residence_energy_per_m2 <- 1.0;
	bool invert_kendall_income <- true;
	float distance_from_object_block_to_grid_km <- 5.0;
	// less commuting criteria
	float min_commute_distance_before_criteria <- 2.5;
	float min_decrease_ratio_criteria <- 0.3;
	float min_decrease_distance_criteria <- 0.5;

	
	// preference parameters
	float b_move_low_inc <- -1.3;
	float b_move_low_inc_kendall <- -1.3;
	float b_commute_distance_low_inc <- -0.005; //-1.28 /-0.88/-0.01
	float b_large_size_low_inc <- 0.0; //0.24
	float b_price_low_inc <- -0.003;    //-0.004
	float b_pop_density_low_inc <- -0.000; //-0.0056
	float b_inc_disparity_low_inc <- -0.0;  //-0.21
	
	float b_move_high_inc <- -1.3;  //-1.29
	float b_move_high_inc_kendall <- -1.6;  //-1.4
	float b_commute_distance_high_inc <- -0.005;  //-1.31 /-0.81/-0.01
	float b_large_size_high_inc <- 0.0;  //0.52
	float b_price_high_inc <- -0.0006;    //-0.0012
	float b_pop_density_high_inc <- -0.0; //-0.0123
	float b_inc_disparity_high_inc <- -0.0;  //-0.52
	
	//policy parameters
	bool incentive_policy <- false;
	bool dynamic_policy <- false;
	float diversity_target <- 0.69;
	float low_inc_pop_ratio_target <- 0.5;
	float commute_distance_decrease_target <- 0.3;
	float building_energy_target <- 55.0;
	float diversity_bottom <- 0.62;
	float low_inc_pop_ratio_bottom <- 0.35;
	float commute_distance_decrease_bottom <- -0.1;
	float building_energy_bottom <- 60.0;
	float max_rent_discount_ratio <- 0.7;
	string new_construction_grids_string <- 'none';
	string new_construction_far_string <- '*1, +0';
	string new_construction_scale <- 'Small';
	float construction_intensity <- 1.0;
	string rent_discount_grids_string <- 'none';
	float rent_discount_ratio_all <- 1.0;
	float rent_discount_ratio_low_inc <- 1.0;
	float rent_discount_ratio_less_commuting <- 1.0;
	float rent_discount_ratio_small_scale <- 1.0;
	float min_finance_to_start_new_construction <- 0.0; //40000.0
	float gap_low_inc_pop_ratio <- 0.0;
	float gap_diversity <- 0.0;
	float gap_commute_distance <-0.0;
	float gap_building_energy <- 0.0;
	list<landuse> incentive_policy_grids <- [];

	// performance indices
	map<string, float> kendall_occupancy <- ['overall'::0.0 ,'small'::0, 'large'::0];
	float kendall_diversity;
	float kendall_low_inc_ratio;
	float kendall_building_energy;
	float kendall_density;
	float kendall_finance;
	map<int,int> move_to_kendall_pop;
	map<int,int> move_out_of_kendall_pop;
	list<int> move_max_value_record_list; //record the max value of move_to_kendall & move_out_of_kendall, so as to set crt_move_y_axis_upper_range
	int crt_move_y_axis_upper_range;   //to set the max value of y-axis range
	map<string, float> mean_commute_distance;
	map<string, float> commute_distance_decrease <- ['totoal'::0.0, 'mean'::0.0];
	map<string, float> kendall_resident_utility;
	map<string, float> kendall_resident_utility_at_start; //use to normalization in viz;
	float residence_energy_per_person;
	
	// viz control
	bool focus_on_kendall <- false;
	bool is_focused <- false;
	bool show_legend <- true;
	bool show_workplace <- false;
	bool show_commute <- false;
	bool show_blockgroup <- true;
	bool show_building <- false;
	bool show_move <- false;
	string show_people <- "Show all people";
	rgb legend_border_color <- #black;
	
	init {
		seed <- 10.0;
		// load some pre-cooked data
		if file_exists(blockgroup_distance_matrix_json_path) {
			float t1 <- machine_time;
			file blockgroup_distance_matrix_json <- json_file(blockgroup_distance_matrix_json_path);
			blockgroup_distance_matrix_map <- blockgroup_distance_matrix_json.contents;
			float t2 <- machine_time;
			write "Load block-to-block distance matrix: " + string((t2-t1)/1000) + " seconds elasped";
		} else {
			write "Prepare block-to-block distance matrix, it may take several minutes";
			do create_blockgroup_dist_matrix;     // takes too many time (5-10 mins) in GAMA...
		}
		
		if file_exists(blockgroup_to_grid_distance_matrix_json_path) {
			float t1 <- machine_time;
			file blockgroup_to_grid_distance_matrix_json <- json_file(blockgroup_to_grid_distance_matrix_json_path);
			blockgroup_to_grid_distance_matrix_map <- blockgroup_to_grid_distance_matrix_json.contents;
			float t2 <- machine_time;
			write "Load block-to-grid distance matrix: " + string((t2-t1)/1000) + " seconds elasped";
		} else {
			write "Prepare block-to-grid distance matrix, it may take several minutes";
			do create_blockgroup_to_grid_dist_matrix;
		}
		
		if file_exists(grid_distance_matrix_json_path) {
			float t1 <- machine_time;
			file grid_distance_matrix_json <- json_file(grid_distance_matrix_json_path);
			grid_distance_matrix_map <- grid_distance_matrix_json.contents;
			float t2 <- machine_time;
			write "Load grid-to-grid distance matrix: " + string((t2-t1)/1000) + " seconds elasped";
		} else {
			write "Prepare grid-to-grid distance matrix, it may take several minutes";
			do create_grid_dist_matrix;
		}
		
		loop gridid over: grid_distance_matrix_map.keys {
			grid_neighbors_map[gridid] <- grid_distance_matrix_map[gridid].keys where (grid_distance_matrix_map[gridid][each] <= 300);
		}

		// create agents
		create blockgroup from: block_shp with: [
			init_low_inc_pop::int(read('low_inc')), 
			init_high_inc_pop::int(read('high_inc')),
			crt_nb_available_bedrooms::int(read('nb_bedroom')),
//			rent::float(read('rent')),
			rent:: max(1000.0, min(1700.0, float(read('rent')))),   //the original distribution of rent is too wide
			geoid::read('GEOID'),
			is_abstract_whole_kendall::false,
			myarea::float(read('real_area'))
		] {
			blockgroup_lookup_map << geoid::self;
			if init_low_inc_pop + init_high_inc_pop > 0 {
				small_size_apt_ratio <- init_low_inc_pop / (init_low_inc_pop + init_high_inc_pop);
			} else {
				small_size_apt_ratio <- 0.0;
			}
			// for kendall blocks, invert low inc an high inc pop so that currently high inc pop are on a larger proportion
			if geoid in kendall_geoid_list {
				int tmp1 <- init_low_inc_pop;
				init_low_inc_pop <- init_high_inc_pop;
				init_high_inc_pop <- tmp1;
			}
		}
		
		create landuse from: grid_shp with: [
			id:: 'g'+read('id'),
			usage:: read('landuse'),
			far:: float(read('FAR')) * far_discount,
			max_height:: int(read('max_height')),
			scale:: read('scale'),
			base_area:: float(read('build_area')),
//			rent_base:: float(read('rent')),
			rent_base:: max(1000.0, min(1700.0, float(read('rent')))),    //the original distribution of rent is too wide
			init_vacant_nb_bedrooms:: int(read('nb_bedroom'))   // not based on landuse(residence capacity) and population, but from appartment shapefile
		] {
			if usage = 'R' and init_vacant_nb_bedrooms = 0 { // do not have rent sample info from appartment shapefile, use mean rent
				rent_base <- scale='S' ? mean_monthly_rent_low_income : mean_monthly_rent_high_income;
			}
			do get_associated_blockgroup;
			do set_rent_policy ();   // update rent_subgroups with all rent_discount_ratio=1 and all rent_shift=0
		}
		kendall_envelope <- envelope(landuse) * 1;
		
//		create building from:buildings_shp;
		create government number:1;
		create developer number:1 returns:tmp;
		the_developer <- first(tmp);
		
		ask blockgroup where (population_discount contains_key each.geoid) {
			init_low_inc_pop <- int(init_low_inc_pop * population_discount[geoid]);
			init_high_inc_pop <- int(init_high_inc_pop * population_discount[geoid]);
		}
		kendall_blockgroups <- blockgroup where (each.geoid in kendall_geoid_list);
		nonkendall_blockgroups <- list(blockgroup) - kendall_blockgroups;
		nonkendall_geoid_list <- nonkendall_blockgroups collect each.geoid;
		
		do load_worker_data;
		
		
		
		do create_people;
//		do set_workplace;   //already embeded into create_people
		
		ask landuse {do update_current_population;}
		ask blockgroup {do update_current_population;}
		
		create blockgroup number:1 returns: kendall_virual_block_list{    //create a virtual blockgroup to represent the whole kendall (all the landuse grids)
			is_abstract_whole_kendall <- true;
			geoid <- 'kendall';
			blockgroup_lookup_map << geoid::self;
			/*
			// use landuse information instead
			 myarea <- sum(kendall_blockgroups collect each.myarea);
			attached_residents <- kendall_blockgroups accumulate each.attached_residents;
			crt_nb_available_bedrooms <- sum(kendall_blockgroups collect each.crt_nb_available_bedrooms);
			if crt_nb_available_bedrooms > 0 {
				rent <- sum(kendall_blockgroups collect (each.rent * each.crt_nb_available_bedrooms)) / crt_nb_available_bedrooms;
			} else {
				rent <- 0.0;
			}
			 */
			myarea <- sum(landuse collect each.myarea);
			attached_residents <- landuse accumulate each.attached_residents;
			crt_nb_available_bedrooms <- sum(landuse collect each.crt_nb_available_bedrooms);
			do get_rent_subgroups_for_kendall;
			rent_subgroups <- ['all_others'::rent, 'low_income'::rent,'less_commuting'::rent, 'low_income_and_less_commuting'::rent];
			nb_workers_in_kendall['low'] <- sum(kendall_blockgroups collect each.nb_workers_in_kendall['low']);
			nb_workers_in_kendall['high'] <- sum(kendall_blockgroups collect each.nb_workers_in_kendall['high']);
			init_low_inc_pop <- sum(kendall_blockgroups collect each.init_low_inc_pop);
			init_high_inc_pop <- sum(kendall_blockgroups collect each.init_high_inc_pop);
			int crt_nb_available_bedrooms_small <- sum((landuse where (each.is_affordable)) collect each.crt_nb_available_bedrooms);
			int crt_nb_available_bedrooms_large <- sum((landuse where (each.is_affordable=false and each.usage='R')) collect each.crt_nb_available_bedrooms);
			write string(crt_nb_available_bedrooms) + ', S=' + string(crt_nb_available_bedrooms_small) + ', L=' + string(crt_nb_available_bedrooms_large);
			if crt_nb_available_bedrooms > 0 {
				small_size_apt_ratio <- crt_nb_available_bedrooms_small / crt_nb_available_bedrooms;
			} else {
				small_size_apt_ratio <- 0.0;
			}
		}
		kendall_virtual_block <- first(kendall_virual_block_list);
		ask landuse {do update_current_population;}
		ask blockgroup {do update_current_population;}
		the_developer.total_rent_at_start <- sum((people where each.live_in_kendall) collect (each.myrent * each.represented_nb_people))*start_rent_coeff;
		do kendall_statistics;
	}
	
	
	reflex main {
		write "\nNew step start:";
		move_to_kendall_pop <- map([0::0, 1::0]);
		move_out_of_kendall_pop <- map([0::0, 1::0]);
		move_stat <- [];
		
		// apply static incentive policy  
		if incentive_policy=true and dynamic_policy=false {
			if cycle = 1 {do apply_static_incentive_policy;}
		}
		
		// apply dynamic incentive policy
		if incentive_policy=true and dynamic_policy=true and mod(cycle,2) = 1 {
			do apply_dynamic_incentive_policy;
		}

		ask people {
			settled_time <- settled_time + 1;   
		}
		
		list<people> all_people <- people where (each.represented_nb_people>0 and each.settled_time>=5);
		ask all_people {less_commuting_if_move_to_kendall <- false;}
		list<people> people_not_in_kendall_and_has_potential_for_less_commuting;
		if (incentive_policy=true and dynamic_policy=false and rent_discount_ratio_less_commuting<0.96) or 
		   (incentive_policy=true and dynamic_policy=true and (gap_commute_distance>0.15)){
			people_not_in_kendall_and_has_potential_for_less_commuting <- people where (
				each.represented_nb_people>0 and each.settled_time>=5 and each.live_in_kendall=false and
				each.commute_distance/1000 >= min_commute_distance_before_criteria and
				each.distance_to_blockgroups_from_my_work_map['250173523001'] < each.commute_distance*(1-min_decrease_ratio_criteria)
			); 
			ask people_not_in_kendall_and_has_potential_for_less_commuting {less_commuting_if_move_to_kendall <- true;}
//			write "total number of less commuting people: " + sum(people_not_in_kendall_and_has_potential_for_less_commuting collect each.represented_nb_people);
			all_people <- people_not_in_kendall_and_has_potential_for_less_commuting + (
				all_people - people_not_in_kendall_and_has_potential_for_less_commuting
			);
		}
//		write "len="+length(people_not_in_kendall_and_has_potential_for_less_commuting);
		ask all_people {
//			write "commute_distance: "+commute_distance + ", to_kendall: " + distance_to_blockgroups_from_my_work_map['250173523001'];
//			if ((work_in_kendall=true and flip(0.6)) or flip(ratio_of_people_considering_move_at_each_step)) {
			float tmp_p <- less_commuting_if_move_to_kendall ? 
				(ratio_of_people_considering_move_at_each_step/(rent_discount_ratio_less_commuting)^3) : ratio_of_people_considering_move_at_each_step;
			if flip(tmp_p) {
//			if flip(1) {
				do select_residence;
				do update_landuse_and_blockgroup_population_info (true, true);
			}
		}
//		write "Structure of move-in people in terms of rent type: ";
//		do move_in_type_statistics;
		ask landuse {
			do update_current_population;
			if usage='R' {do update_attached_resident_rent;}
		}
		ask blockgroup {
			do update_current_population;
		}
//		write 'update over, vacant house=' + kendall_virtual_block.crt_nb_available_bedrooms;
		do kendall_statistics;
		ask the_developer {do collect_additional_rent_and_pay_subsidy;}
		ask people where (each.represented_nb_people <= 0) {
			if represented_nb_people > 0 {
				show_size <- 5 + (20-5)*(represented_nb_people-1)/(50-1);
			} else {
				show_size <- 0.0;
			}
			
		}
		do reference_info;
	}
	
	reflex do_focus when: is_focused != focus_on_kendall{
		if focus_on_kendall{
			focus_on kendall_envelope;
			is_focused <- true;
			show_people <- 'Show Kendall people only';
		} else {
			focus_on shape;
			is_focused <- false;
			show_people <- 'Show all people';
		}
	}
	
	reflex halting when: cycle>max_cycle {
		do pause;
	}
	
	action update_landuse_and_blockgroup_population_info (bool only_kendall <- false, bool quick <- false) {
//		float t1 <- machine_time;
		if quick {
			ask landuse where (each.usage='R') {do update_current_population (true);}
			if only_kendall {
				ask kendall_virtual_block {do update_current_population (true);}
			} else {
				ask blockgroup {do update_current_population (true);}
			}
		} else {
			ask landuse {do update_current_population (false);}
			if only_kendall {
				ask kendall_virtual_block {do update_current_population (false);}
			} else {
				ask blockgroup {do update_current_population (false);}
			}
		}
//		float t2 <- machine_time;
//		write 'update population info cost time: ' + (t2-t1);
	}
	
	action reference_info {
//		write machine_time;
//		list<people> tmp_people <- people where (each.work_in_kendall and each.live_in_kendall=false and each.represented_nb_people>0);
		
//		write "Debug info start: \n";
//		write '\n\n==================\nPoeple work in Kendall but live outside: ' + sum(tmp_people collect each.represented_nb_people);
//		write tmp_people;
//		write tmp_people collect each.represented_nb_people;

//		int nb_all_others;
//		int nb_less_commuting;
//		int nb_less_commuting_and_low_income;
//		int nb_low_income;
//		int nb_live_in_kendall;
//		ask people where each.live_in_kendall {
//			nb_live_in_kendall <- nb_live_in_kendall + represented_nb_people;
//			if my_rent_type = 'all_others' {nb_all_others <- nb_all_others + represented_nb_people;}
//			else if my_rent_type = 'less_commuting' {nb_less_commuting <- nb_less_commuting + represented_nb_people;}
//			else if my_rent_type = 'low_income_and_less_commuting' {nb_less_commuting_and_low_income <- nb_less_commuting_and_low_income + represented_nb_people;}
//			else if my_rent_type = 'low_income' {nb_low_income <- nb_low_income + represented_nb_people;}
//			else {write "Unknow rent_type: " + my_rent_type;}
//		}
//		write  'total=' + nb_live_in_kendall + ', all_others=' + nb_all_others + ', low_income=' + nb_low_income + ', less_commuting=' + 
//			nb_less_commuting + ', low_income_and_less_commuting=' + nb_less_commuting_and_low_income;
//		write kendall_virtual_block.rent_subgroups;
//		write 'mean_commute_distance: '+ mean_commute_distance['live_in_kendall'];
//		write nb_nonkendall_blocks_in_hom_loc_choice/length(nonkendall_geoid_list+1);
//		write kendall_virtual_block.crt_total_pop;
//		write "kendall_diversity: "+kendall_diversity+", kendall_low_inc_ratio: "+kendall_low_inc_ratio;
	}
	
	action load_worker_data {
		matrix kendall_workers <- matrix(workers_csv);
		loop idx from:0 to: kendall_workers.rows-1 {
			string this_geoid <- kendall_workers[0, idx];
			blockgroup residence_blockgroup <- first(blockgroup where (each.geoid=this_geoid));
			if residence_blockgroup=nil {residence_blockgroup <- any(kendall_blockgroups);}
			int this_nb_works <- int(kendall_workers[1, idx]);
			string this_income <- kendall_workers[2, idx];
			residence_blockgroup.nb_workers_in_kendall[this_income] <- residence_blockgroup.nb_workers_in_kendall[this_income] + this_nb_works;
		}
	}
	
	action create_people {
		list<blockgroup> all_blockgroups <- list(blockgroup);
		ask blockgroup where (each.is_abstract_whole_kendall = false) {
			int this_represented_nb_people;
			bool this_in_kendall <- (geoid in kendall_geoid_list) ? true : false;
			list<people> this_low_income_population_work_in_kendall <- [];
			list<people> this_high_income_population_work_in_kendall <- [];
			list<people> this_low_income_population_work_outof_kendall <- [];
			list<people> this_high_income_population_work_outof_kendall <- [];
			
			this_represented_nb_people <- this_in_kendall ? lower_nb_represented_persons : upper_nb_represented_persons;
			list<people> this_low_income_population <- create_people_given_income_and_blockgroup  (
				int((init_low_inc_pop) / this_represented_nb_people), 0, 
				this_in_kendall, this_represented_nb_people, #red
			);
			list<people> this_high_income_population <- create_people_given_income_and_blockgroup  (
				int((init_high_inc_pop) / this_represented_nb_people), 1, 
				this_in_kendall, this_represented_nb_people, #blue
			);
			ask this_low_income_population + this_high_income_population {
				work_in_kendall <- false;
				settled_time <- 10000;
				// randomly determin workplace using a very simple logit model with commute distance as the single regressor
				list<float> dist_to_blockgruops_from_my_home <- all_blockgroups collect (blockgroup_distance_matrix_map[home_block.geoid][each.geoid] / 1000);
				list<float> exp_v;
				if this_in_kendall {
					exp_v <- dist_to_blockgruops_from_my_home collect exp(each* beta_work_allocation_for_kendall_residents);
				} else {
					exp_v <- dist_to_blockgruops_from_my_home collect exp(each* beta_work_allocation);
				}
				int choice_idx <- rnd_choice(exp_v);
				work_block <- all_blockgroups[choice_idx];
				if work_block.geoid in kendall_geoid_list {
					work_in_kendall <- true;
					list<landuse> work_grid_candidates <- landuse where (each.worker_capacity - each.crt_nb_workers >= represented_nb_people);
					if length(work_grid_candidates) > 0 {
						work_grid <- any(work_grid_candidates);
						work_grid.crt_nb_workers <- work_grid.crt_nb_workers + represented_nb_people;
						work_grid.attached_workers << self;
						work_loc <- any_location_in(work_grid);
						work_block <- work_grid.associated_blockgroup;
					} else {
						write "Set workplace for Kendall workers: cannot find a inside workplace";
						work_block <- any(blockgroup);
						work_loc <- any_location_in(work_block);
					}
				} else {
					work_loc <- any_location_in(work_block);
				}
				work_block.attached_workers << self;
				
				distance_to_blockgroups_from_my_work_map <- blockgroup_distance_matrix_map[work_block.geoid];
				nonkendall_geoids_in_3km_from_work <- nonkendall_geoid_list where (distance_to_blockgroups_from_my_work_map[each] <= 3000.0);
				if work_grid != nil {
					distance_to_grids_from_my_work_map <- grid_distance_matrix_map[work_grid.id];
				} else {
					distance_to_grids_from_my_work_map <- blockgroup_to_grid_distance_matrix_map[work_block.geoid];	
				}
				if work_grid != nil and home_grid != nil {commute_distance <- grid_distance_matrix_map[home_grid.id][work_grid.id];}
				else if work_grid = nil and home_grid = nil {commute_distance <- blockgroup_distance_matrix_map[home_block.geoid][work_block.geoid];}
				else if work_grid != nil {commute_distance <- blockgroup_to_grid_distance_matrix_map[home_block.geoid][work_grid.id];}
				else if home_grid != nil {commute_distance <- blockgroup_to_grid_distance_matrix_map[work_block.geoid][home_grid.id];}
				commute_distance_init <- commute_distance;
			}
			
			// ...
			/* 
			// create people who work in Kendall
			if nb_workers_in_kendall['low'] > 0 {
				this_represented_nb_people <- nb_workers_in_kendall['low'] >= lower_nb_represented_persons ? lower_nb_represented_persons : 1;
				this_low_income_population_work_in_kendall <- create_people_given_income_and_blockgroup (
					nb_workers_in_kendall['low'], 0, this_in_kendall, this_represented_nb_people, #red
				);
			} 
			if nb_workers_in_kendall['high'] > 0 {
				this_represented_nb_people <- nb_workers_in_kendall['high'] >= lower_nb_represented_persons ? lower_nb_represented_persons : 1;
				this_high_income_population_work_in_kendall <- create_people_given_income_and_blockgroup (
					nb_workers_in_kendall['high'], 1, this_in_kendall, this_represented_nb_people, #blue
				);
			}
			ask this_low_income_population_work_in_kendall + this_high_income_population_work_in_kendall {
				work_in_kendall <- true;
				settled_time <- 10000;
				list<landuse> work_grid_candidates <- landuse where (each.worker_capacity - each.crt_nb_workers >= represented_nb_people);
				if length(work_grid_candidates) > 0 {
					work_grid <- any(work_grid_candidates);
					work_grid.crt_nb_workers <- work_grid.crt_nb_workers + represented_nb_people;
					work_grid.attached_workers << self;
					work_loc <- any_location_in(work_grid);
					work_block <- work_grid.associated_blockgroup;
				} else {
					write "Set workplace for Kendall workers: cannot find a inside workplace";
					work_block <- any(blockgroup);
					work_loc <- any_location_in(work_block);
				}
				distance_to_blockgroups_from_my_work_map <- blockgroup_distance_matrix_map[work_block.geoid];
				nonkendall_geoids_in_3km_from_work <- nonkendall_geoid_list where (distance_to_blockgroups_from_my_work_map[each] <= 3000.0);
				if work_grid != nil {
					distance_to_grids_from_my_work_map <- grid_distance_matrix_map[work_grid.id];
				} else {
					distance_to_grids_from_my_work_map <- blockgroup_to_grid_distance_matrix_map[work_block.geoid];	
				}
				if work_grid != nil and home_grid != nil {commute_distance <- grid_distance_matrix_map[home_grid.id][work_grid.id];}
				else if work_grid = nil and home_grid = nil {commute_distance <- blockgroup_distance_matrix_map[home_block.geoid][work_block.geoid];}
				else if work_grid != nil {commute_distance <- blockgroup_to_grid_distance_matrix_map[home_block.geoid][work_grid.id];}
				else if home_grid != nil {commute_distance <- blockgroup_to_grid_distance_matrix_map[work_block.geoid][home_grid.id];}
			}
			
			// create people who work out of Kendall
			this_represented_nb_people <- this_in_kendall ? lower_nb_represented_persons : upper_nb_represented_persons;
			this_low_income_population_work_outof_kendall <- create_people_given_income_and_blockgroup  (
				int((init_low_inc_pop - nb_workers_in_kendall['low']) / this_represented_nb_people), 0, 
				this_in_kendall, this_represented_nb_people, #red
			);
			this_high_income_population_work_outof_kendall <- create_people_given_income_and_blockgroup  (
				int((init_high_inc_pop - nb_workers_in_kendall['high']) / this_represented_nb_people), 1, 
				this_in_kendall, this_represented_nb_people, #blue
			);
			ask this_low_income_population_work_outof_kendall + this_high_income_population_work_outof_kendall {
				work_in_kendall <- false;
				settled_time <- 10000;
				// randomly determin workplace using a very simple logit model with commute distance as the single regressor
				list<float> dist_to_nonkendall_blockgruops_from_my_home <- nonkendall_blockgroups collect (blockgroup_distance_matrix_map[home_block.geoid][each.geoid] / 1000);
				list<float> exp_v;
				if this_in_kendall {
					exp_v <- dist_to_nonkendall_blockgruops_from_my_home collect exp(each* beta_work_allocation_for_kendall_residents);
				} else {
					exp_v <- dist_to_nonkendall_blockgruops_from_my_home collect exp(each* beta_work_allocation);
				}
				int choice_idx <- rnd_choice(exp_v);
				work_block <- nonkendall_blockgroups[choice_idx];
//				work_block <- any(nonkendall_blockgroups);
				work_loc <- any_location_in(work_block);
				distance_to_blockgroups_from_my_work_map <- blockgroup_distance_matrix_map[work_block.geoid];
				nonkendall_geoids_in_3km_from_work <- nonkendall_geoid_list where (distance_to_blockgroups_from_my_work_map[each] <= 3000.0);
				if work_grid != nil {
					distance_to_grids_from_my_work_map <- grid_distance_matrix_map[work_grid.id];
				} else {
					distance_to_grids_from_my_work_map <- blockgroup_to_grid_distance_matrix_map[work_block.geoid];	
				}
				if work_grid != nil and home_grid != nil {commute_distance <- grid_distance_matrix_map[home_grid.id][work_grid.id];}
				else if work_grid = nil and home_grid = nil {commute_distance <- blockgroup_distance_matrix_map[home_block.geoid][work_block.geoid];}
				else if work_grid != nil {commute_distance <- blockgroup_to_grid_distance_matrix_map[home_block.geoid][work_grid.id];}
				else if home_grid != nil {commute_distance <- blockgroup_to_grid_distance_matrix_map[work_block.geoid][home_grid.id];}
			}
			*/
		}
	}
	
	action set_workplace {
		// work in Kendall
		matrix kendall_workers <- matrix(workers_csv);
		loop idx from:0 to: kendall_workers.rows-1 {
			string this_geoid <- kendall_workers[0, idx];
			blockgroup residence_blockgroup <- first(blockgroup where (each.geoid=this_geoid));
			if residence_blockgroup = nil {residence_blockgroup <- any(blockgroup);}
			int this_nb_works <- int(kendall_workers[1, idx]);
			int this_income_type <- kendall_workers[2, idx] = "low" ? 0 : 1;
			rgb this_color <- this_income_type=0 ? #red : #blue;
			bool this_in_kendall <- this_geoid in kendall_geoid_list;
			float this_show_size <- 5 + (20-5)*(1-1)/(50-1);
			
			ask this_nb_works among (residence_blockgroup.attached_residents) {
				list<landuse> work_grid_candidates <- landuse where (each.worker_capacity >= represented_nb_people);
				if length(work_grid_candidates) > 0 {
					work_grid <- any(work_grid_candidates);
					work_loc <- any_location_in(work_grid);
					work_block <- work_grid.associated_blockgroup;
				} else {
					write "Set workplace for Kendall workers: cannot find a inside workplace";
					work_block <- any(blockgroup);
					work_loc <- any_location_in(work_block);
				}
			}
		}
		// work ouside Kendall
		ask people where (each.work_loc=nil) {
			work_block <- any(nonkendall_blockgroups);
			work_loc <- any_location_in(work_block);
		}
	}
	
	
	action kendall_statistics {
		// occupany
		kendall_occupancy['overall'] <- kendall_virtual_block.crt_total_pop / sum(landuse collect each.resident_capacity);
		list<landuse> affordable_residence_grids <- landuse where (each.usage='R' and each.is_affordable);
		list<people> people_in_affordable_residences <- affordable_residence_grids accumulate each.attached_residents;
		list<landuse> normal_residence_grids <- landuse where (each.usage='R' and each.is_affordable=false);
		list<people> people_in_normal_residences <- normal_residence_grids accumulate each.attached_residents;
		kendall_occupancy['small'] <- sum(people_in_affordable_residences collect each.represented_nb_people) / 
			sum(affordable_residence_grids collect each.resident_capacity);
		kendall_occupancy['large'] <- sum(people_in_normal_residences collect each.represented_nb_people) / 
			sum(normal_residence_grids collect each.resident_capacity);
		
		kendall_diversity <- kendall_virtual_block.crt_diversity;
		kendall_low_inc_ratio <- kendall_virtual_block.crt_low_inc_pop / kendall_virtual_block.crt_total_pop;
		kendall_building_energy <- 0.0;
		kendall_density <- kendall_virtual_block.crt_pop_density;
		kendall_finance <- 0.0;
		
		// commute calculation
//		mean_commute_distance['all'] <- sum(people collect (each.represented_nb_people * each.commute_distance)) / sum(people collect each.represented_nb_people);
		list<people> people_all <- people where (each.represented_nb_people>0);
		list<people> people_work_in_kendall <- people where (each.work_in_kendall and each.represented_nb_people>0);
		list<people> people_live_in_kendall <- people where (each.live_in_kendall and each.represented_nb_people>0);
		list<people> people_work_or_live_in_kendall <- people where ((each.work_in_kendall or each.live_in_kendall) and each.represented_nb_people>0);
		mean_commute_distance['all'] <- sum(people_all collect (each.represented_nb_people * each.commute_distance)) / 
			sum(people_all collect each.represented_nb_people);
		mean_commute_distance['work_in_kendall'] <- sum(people_work_in_kendall collect (each.represented_nb_people * each.commute_distance)) / 
			sum(people_work_in_kendall collect each.represented_nb_people);
		mean_commute_distance['live_in_kendall'] <- sum(people_live_in_kendall collect (each.represented_nb_people * each.commute_distance)) / 
			sum(people_live_in_kendall collect each.represented_nb_people);
		mean_commute_distance['work_or_live_in_kendall'] <- sum(people_work_or_live_in_kendall collect (each.represented_nb_people * each.commute_distance)) / 
			sum(people_work_or_live_in_kendall collect each.represented_nb_people);
			
		float total_commute_distance_decrease <- sum((people where each.live_in_kendall) collect (
			each.represented_nb_people * (each.commute_distance_init - each.commute_distance) / 1000
		));
		float mean_commute_distance_decrease <- total_commute_distance_decrease / sum((people where each.live_in_kendall) collect each.represented_nb_people);
		commute_distance_decrease <- ['total'::total_commute_distance_decrease, 'mean'::mean_commute_distance_decrease];
			
		move_max_value_record_list << max(move_to_kendall_pop.values + move_out_of_kendall_pop);  
//		move_max_value_record_list <- [1,2,3,4,5,6,7,8] collect move_max_value_record_list[max(length(move_max_value_record_list)-each, 0)];  //keep only the last 3 elements
//		crt_move_y_axis_upper_range <- int(quantile(move_max_value_record_list, 0.999999999999999))+10;
		if cycle=4 {remove max(move_max_value_record_list) from:move_max_value_record_list;}
		crt_move_y_axis_upper_range <- max(move_max_value_record_list) + 15;
//		write string(crt_move_y_axis_upper_range) + ', ' +  string(int(quantile(move_max_value_record_list, 1))+10) + ', ' + string(max(move_max_value_record_list));
		
		// utility calculation
		ask people where (each.live_in_kendall) {
			do calculate_utility;
		}
		kendall_resident_utility['total'] <- sum( (people where (each.live_in_kendall)) collect (each.represented_nb_people * each.utility) ) / 
			sum((people where (each.live_in_kendall)) collect each.represented_nb_people);
		kendall_resident_utility['low_inc'] <- sum( (people where (each.live_in_kendall and each.income=0)) collect (each.represented_nb_people * each.utility) ) / 
			sum((people where (each.live_in_kendall and each.income=0)) collect each.represented_nb_people);
		kendall_resident_utility['high_inc'] <- sum( (people where (each.live_in_kendall and each.income=1)) collect (each.represented_nb_people * each.utility) ) / 
			sum((people where (each.live_in_kendall and each.income=1)) collect each.represented_nb_people);	
		if cycle = 0 {
			kendall_resident_utility_at_start <- ['total'::kendall_resident_utility['total'], 'low_inc'::kendall_resident_utility['low_inc'], 'high_inc'::kendall_resident_utility['high_inc']];
		}
//		write kendall_resident_utility;

		//energy
		float total_occupied_area_small_scale <- sum((landuse where (each.scale='S')) collect each.crt_total_pop) * residence_area_small_scale;
		float total_occupied_area_large_scale <- sum((landuse where (each.scale!='S')) collect each.crt_total_pop) * residence_area_large_scale;
		residence_energy_per_person <-  residence_energy_per_m2 * (total_occupied_area_small_scale + total_occupied_area_large_scale) / kendall_virtual_block.crt_total_pop;
	}
	
	map<string, int> move_in_type_statistics (int max_settled_time<-0, bool printFlag<-true) {
		map<string,int> results;
		list<people> move_in_people <- people where (each.live_in_kendall and each.represented_nb_people>0 and each.settled_time<=max_settled_time);
		loop this_rent_type over: ['low_income', 'less_commuting', 'low_income_and_less_commuting', 'all_others'] {
			results[this_rent_type] <- sum((move_in_people where (each.my_rent_type=this_rent_type)) collect each.represented_nb_people);
		}
		if sum(move_in_people collect each.represented_nb_people) > sum(results.values) {
			results['error'] <- sum(move_in_people collect each.represented_nb_people) - sum(results.values);
		}
		if printFlag {
			write results;
		}
		return results;
	}
	
	action apply_static_incentive_policy {
//		list<landuse> rent_discount_grids <- 5 among (landuse where (each.usage='R'));
		list<landuse> rent_discount_grids <- 6 among (landuse where (each.usage='vacant'));
//		list<landuse> rent_discount_grids <- landuse where (each.id in ['g80', 'g79', 'g81', 'g82', 'g83']);
//		list<landuse> new_construction_grids <- 5 among (landuse where ((each.usage='R' and each.scale='R') or each.usage='vacant'));
		list<landuse> new_construction_grids <- rent_discount_grids;
		//calculate subsidy need to pay due to rent discount
		ask rent_discount_grids {
			do set_rent_policy (['all'::rent_discount_ratio_all, 'low_income'::rent_discount_ratio_low_inc,
				 'small_scale'::rent_discount_ratio_small_scale, 'less_commuting'::rent_discount_ratio_less_commuting]
			);
		}
//		float expected_subsidy <- rent_discount_grids collect (each.crt_nb_available_bedrooms*0.8*each.rent_base);
		float expected_subsidy <- 0.0;
		ask rent_discount_grids {
			expected_subsidy <- expected_subsidy + crt_nb_available_bedrooms * (rent_base - min(rent_subgroups.values)) * 0.7;
		}
		
		//anyway, do it first
//		float new_construction_far_shift_max <- 1.2;
		float total_rent_discount_ratio <- rent_discount_ratio_low_inc * rent_discount_ratio_small_scale * rent_discount_ratio_less_commuting;
		float new_construction_far_shift_max <- new_construction_far_given_incentive_policy(total_rent_discount_ratio);
		float new_construction_far_shift <- new_construction_far_shift_max * construction_intensity;
		ask new_construction_grids {
			do new_constructions(1.0, new_construction_far_shift);
//			do change_landuse('R', 'S');
		}
		int nb_small_scale_for_new_construction <- 0;
		if rent_discount_ratio_small_scale <= 0.5 {
			nb_small_scale_for_new_construction <- 6;
		} else if rent_discount_ratio_small_scale <= 0.7 {
			nb_small_scale_for_new_construction <- 4;
		} else if rent_discount_ratio_small_scale <= 0.95 {
			nb_small_scale_for_new_construction <- 2;
		} else {
			nb_small_scale_for_new_construction <- 0;
		}
		list<landuse> new_construction_grids_small_scale <- nb_small_scale_for_new_construction among new_construction_grids;
		list<landuse> new_construction_grids_large_scale <- new_construction_grids - new_construction_grids_small_scale;
		float rent_for_new_constructions <- 1100.0;
		ask new_construction_grids_small_scale {
			do change_landuse('R', 'S');
			rent_base <- rent_for_new_constructions;
		}
		ask new_construction_grids_large_scale {
			do change_landuse('R', 'L');
			rent_base <- rent_for_new_constructions;
		}
		
		ask rent_discount_grids {
			do set_rent_policy (['all'::rent_discount_ratio_all, 'low_income'::rent_discount_ratio_low_inc,
				 'small_scale'::rent_discount_ratio_small_scale, 'less_commuting'::rent_discount_ratio_less_commuting]
			);
		}
		
		
		
	}
	
	action apply_static_incentive_policy_complex{
		list<landuse> new_construction_grids <- parse_input_grids(new_construction_grids_string);
		list<landuse> rent_discount_grids <- parse_input_grids(rent_discount_grids_string);
		list<string> tmp <- new_construction_far_string split_with ',';
		float new_construction_far_multiplier <- float(last(tmp[0] split_with '*'));
		float new_construction_far_shift <- float(last(tmp[1] split_with '+'));
		new_construction_scale <- new_construction_scale='Small' ? 'S' : 'L';
		ask new_construction_grids {
			do new_constructions(new_construction_far_multiplier, new_construction_far_shift);
			do change_landuse('R', new_construction_scale);
		}
		ask rent_discount_grids {
			do set_rent_policy (['all'::rent_discount_ratio_all, 'low_income'::rent_discount_ratio_low_inc,
				 'small_scale'::rent_discount_ratio_small_scale, 'less_commuting'::rent_discount_ratio_less_commuting]
			);
		}
	}
	
	list<landuse> parse_input_grids (string grids_string) {
		grids_string <- lower_case(grids_string);
		list<landuse> result_grids;
		if grids_string = 'none' {
			result_grids <- [];
		} else if grids_string = 'all' {
			result_grids <- landuse where (each.usage in ['vacant', 'R']);
		} else if grids_string = 'vacant' {
			result_grids <- landuse where (each.usage = 'vacant');
		} else if grids_string = 'small' {
			result_grids <- landuse where (each.usage='R' and each.scale='S');
		} else if grids_string = 'large' {
			result_grids <- landuse where (each.usage='R' and each.scale!='S');
		} else if copy_between(grids_string, 0, 3) = 'all' and length(grids_string)>3 {
			int tmp_nb <- int(last(grids_string split_with 'all'));
			result_grids <- tmp_nb among (landuse where (each.usage in ['vacant', 'R']));
		} else if copy_between(grids_string, 0, 6) = 'vacant' and length(grids_string)>6 {
			int tmp_nb <- int(last(grids_string split_with 'vacant'));
			result_grids <- tmp_nb among (landuse where (each.usage='vacant'));
		} else if copy_between(grids_string, 0, 5) = 'small' and length(grids_string)>5 {
			int tmp_nb <- int(last(grids_string split_with 'small'));
			result_grids <- tmp_nb among (landuse where (each.usage='R' and each.scale='S'));
		} else if copy_between(grids_string, 0, 5) = 'large' and length(grids_string)>5 {
			int tmp_nb <- int(last(grids_string split_with 'large'));
			result_grids <- tmp_nb among (landuse where (each.usage='R' and each.scale!='S'));
		} else {
			try {
				list<string> tmp_gridids <- (grids_string split_with ',') collect ('g'+int(each));
				result_grids <- landuse where (each.id in tmp_gridids);
			} catch {
				result_grids <- [];
				write "Cannot parse grids string: " + grids_string;
			}
		}
		return result_grids;
	}
	
	action apply_dynamic_incentive_policy {
		// clear all rent policy before
		ask landuse {
			do set_rent_policy (['low_income':: 1, 'small_scale'::1, 'less_commuting'::1]);
		}
		// gap calculation
		gap_low_inc_pop_ratio <- normalized_gap(kendall_low_inc_ratio, low_inc_pop_ratio_target, low_inc_pop_ratio_bottom);
		write "fuck: " + kendall_low_inc_ratio + ', target = ' + low_inc_pop_ratio_target + ', bottom = ' + low_inc_pop_ratio_bottom;
		gap_diversity <- normalized_gap(kendall_diversity, diversity_target, diversity_bottom);
		gap_commute_distance <- normalized_gap(
			commute_distance_decrease['mean'], commute_distance_decrease_target, commute_distance_decrease_bottom
		);
		gap_building_energy <- normalized_gap(residence_energy_per_person, building_energy_target, building_energy_bottom, false);
		write 'bulding: ' + [residence_energy_per_person, building_energy_target, building_energy_bottom];
		write "Curret gap: \naffordability=" + gap_low_inc_pop_ratio + ', diversity='+gap_diversity + '\ncommute='+gap_commute_distance + ', building='+gap_building_energy;
		ask landuse where (each.usage in ['R', 'vacant']) {
			do target_subgroups_calc(3.0, ['low_income'::max(gap_low_inc_pop_ratio, gap_diversity), 'less_commuting'::gap_commute_distance]);
		}
		// min, max
		list<int> target_low_income_pop_in_all_grids <- landuse where (each.usage in ['R', 'vacant']) collect (each.target_subgroups['low_income']);
		list<int> target_less_commuting_pop_in_all_grids <- landuse where (each.usage in ['R', 'vacant']) collect (each.target_subgroups['less_commuting']);
		int max_target_low_income_pop <- max(target_low_income_pop_in_all_grids);
		int min_target_low_income_po <- min(target_low_income_pop_in_all_grids);
		int max_target_less_commuting_pop <- max(target_less_commuting_pop_in_all_grids);
		int min_target_less_commuting_pop <- min(target_less_commuting_pop_in_all_grids);
//		write "potential calculation\nrent_type_weights = " + ['low_income'::max(gap_low_inc_pop_ratio, gap_diversity), 'less_commuting'::gap_commute_distance];
//		write "low income: max target = " + max_target_low_income_pop + ', min target = ' + min_target_low_income_po;
//		write "less commuting: max target = " + max_target_less_commuting_pop + ', min target = ' + min_target_less_commuting_pop;
		ask landuse where (each.usage in ['R', 'vacant']) {
			do potential_calc(
				['low_income'::max(gap_low_inc_pop_ratio, gap_diversity), 'less_commuting'::gap_commute_distance],
				max_target_low_income_pop, min_target_low_income_po, max_target_less_commuting_pop, min_target_less_commuting_pop
			);
		}
		
//		ask landuse where (each.usage in ['R', 'vacant']) {
//			write 'name = '+name + ', potential = ' + crt_potential + ', target = ' + target_subgroups;
//		}
		list<landuse> grids_sorted_by_potential <- (landuse where (each.usage in ['R', 'vacant'])) sort_by (-each.crt_potential);
		list<landuse> grids_with_top6_potential <- [0,1,2,3,4,5] collect grids_sorted_by_potential[each];
//		list<landuse> grids_sorted_by_target_less_commuting_pop <- (landuse where (each.usage in ['R', 'vacant'])) sort_by (-each.target_subgroups['less_commuting']);
//		list<landuse> grids_with_top6_target_less_commuting_pop <- [0,1,2,3,4,5,6] collect  grids_sorted_by_target_less_commuting_pop[each];
//		list<landuse> grids_sorted_by_low_inc_pop <- (landuse where (each.usage in ['R', 'vacant'])) sort_by (-each.target_subgroups['low_income']);
//		list<landuse> grids_with_top6_target_low_inc_pop <- [0,1,2,3,4,5] collect grids_sorted_by_low_inc_pop[each];
		
		int nb_small_scale_for_new_construction <- 0;
		if gap_building_energy >= 0.7 {
			nb_small_scale_for_new_construction <- 6;
		} else if gap_building_energy >= 0.4 {
			nb_small_scale_for_new_construction <- 6;
		} else if gap_building_energy >= 0.2 {
			nb_small_scale_for_new_construction <- 4;
		} else if gap_building_energy >= 0 {
			nb_small_scale_for_new_construction <- 2;
		} else {
			nb_small_scale_for_new_construction <- 1;
		}
		
		rent_discount_ratio_less_commuting <- 1 - max_rent_discount_ratio * gap_commute_distance;
		rent_discount_ratio_low_inc <- 1 - max_rent_discount_ratio * max(gap_low_inc_pop_ratio, gap_diversity);
		rent_discount_ratio_small_scale <- 1 - max_rent_discount_ratio * gap_building_energy;
		float total_rent_discount_ratio <- rent_discount_ratio_low_inc * rent_discount_ratio_small_scale * rent_discount_ratio_less_commuting;
//		float new_construction_far_shift_max <- new_construction_far_given_incentive_policy(total_rent_discount_ratio);
//		if total_rent_discount_ratio <= 0.4 {
//			float scale_tmp <- (0.4 / (rent_discount_ratio_less_commuting*rent_discount_ratio_low_inc*rent_discount_ratio_small_scale)) ^ (1/3);
//			rent_discount_ratio_less_commuting <- rent_discount_ratio_less_commuting * scale_tmp;
//			rent_discount_ratio_low_inc <- rent_discount_ratio_low_inc * scale_tmp;
//			rent_discount_ratio_small_scale <- rent_discount_ratio_small_scale * scale_tmp;
//		}
		
		
		if the_developer.finance > min_finance_to_start_new_construction and kendall_occupancy['overall']>0.95 {
			ask grids_with_top6_potential {
				do new_constructions(1.0, 1.1);
				if usage='vacant' {
					do change_landuse('R', 'L');
					rent_base <- 1100.0;
				}
			}
			ask nb_small_scale_for_new_construction among grids_with_top6_potential {
				do change_landuse('R', 'S');
				rent_base <- 1100.0;
			}
		}
//		if gap_commute_distance > 0 {
//			ask grids_with_top6_target_less_commuting_pop {
//				do set_rent_policy(['less_commuting':: max_rent_discount_ratio * gap_commute_distance]);
//			}
//		}
//		if gap_low_inc_pop_ratio>0 or gap_diversity>0 {
//			ask grids_with_top6_target_low_inc_pop {
//				do set_rent_policy(['low_income':: max_rent_discount_ratio * max(gap_low_inc_pop_ratio, gap_diversity)]);
//			}
//		}
//		if gap_building_energy > 0 {
//			list<landuse> grids_to_apply_small_scale_discount;
//			if new_construction_this_step = true {
//				grids_to_apply_small_scale_discount <- [];////
//			} else {
//				grids_to_apply_small_scale_discount <- 6 among landuse where (each.usage='R' and each.scale='S');
//			}
//		}
		ask grids_with_top6_potential {
			do set_rent_policy (['low_income':: rent_discount_ratio_low_inc,
				 'small_scale':: rent_discount_ratio_small_scale, 
				 'less_commuting':: rent_discount_ratio_less_commuting]
			); 
		}
		
		write "\n=======================================";
		write "Rent Discount Policy on the following grids: ";
		write grids_with_top6_potential;
		write "Rent Discount for Low Income People: " + rent_discount_ratio_low_inc;
		write "Rent Discount for Small Scale Residence: " + rent_discount_ratio_small_scale;
		write "Rent Discount for Less Commuting People: " + rent_discount_ratio_less_commuting;
		write "=======================================\n";
	}
	
	float normalized_gap (float crt_value, float target_value, float bottom_value, bool max_is_target <- true) {
		float normalized_value;
		if target_value = bottom_value {normalized_value <- 1.0;}
		else if max_is_target = true {
			normalized_value <- (crt_value - bottom_value) / (target_value - bottom_value);
		} else if max_is_target = false {
			normalized_value <- (bottom_value - crt_value) / (bottom_value - target_value);
		}
		return 1 - max(min(normalized_value, 1), 0);
	}
	
	float new_construction_far_given_incentive_policy (float total_rent_discount_ratio) {
		float far_shift <- 0.0;
		float far_shift_allowed <- min(5 * (1 - total_rent_discount_ratio), 1.1);
		float max_debt <- 4000000.0;
		float rent_base_tmp <- 1300.0;
//		float expected_occupancy <- kendall_occupancy['overall'];
		float expected_occupancy <- 1.0;
		float construction_fee_multiplier <- 6.0;
		float constrution_cost_tmp <- construction_fee_multiplier * 240;
		float new_construction_area_accepted <- max_debt / (constrution_cost_tmp - rent_base_tmp/45*3*16 * expected_occupancy * total_rent_discount_ratio);
		float far_shift_accepted <- new_construction_area_accepted / 6400;
		far_shift <- min(far_shift_allowed, far_shift_accepted);
		return far_shift;
	}
	
	action create_blockgroup_dist_matrix {
		float t1 <- machine_time;
		list<string> geoid_list <- blockgroup collect each.geoid;
		loop i_idx from:0 to:(length(geoid_list)-2) {
			string from_geoid <- geoid_list[i_idx];
			blockgroup from_blockgroup <- first(blockgroup where (each.geoid = from_geoid));
			if from_geoid in blockgroup_distance_matrix_map.keys {
				blockgroup_distance_matrix_map[from_geoid] << from_geoid::0.0;
			} else {
				blockgroup_distance_matrix_map[from_geoid] <- map(from_geoid::0.0);
			}
			loop j_idx from:(i_idx+1) to:(length(geoid_list)-1) {
				string to_geoid <- geoid_list[j_idx];
				blockgroup to_blockgroup <- first(blockgroup where (each.geoid = to_geoid));
				float this_dist <- distance_to(centroid(from_blockgroup), centroid(to_blockgroup));
				blockgroup_distance_matrix_map[from_geoid] << to_geoid::this_dist;
				if to_geoid in blockgroup_distance_matrix_map.keys {
					blockgroup_distance_matrix_map[to_geoid] << from_geoid::this_dist;
				} else {
					blockgroup_distance_matrix_map[to_geoid] <- map(from_geoid::this_dist);
				}
			}
		}
		string last_geoid <- geoid_list[length(geoid_list)-1];
		blockgroup_distance_matrix_map[last_geoid] << last_geoid::0;
		file f <- json_file(blockgroup_distance_matrix_json_path, blockgroup_distance_matrix_map);
		save f;
		float t2 <- machine_time;
		write "Create block-block distance matrix: " + string((t2-t1)/1000) + " seconds elasped";
	}
	
	action create_blockgroup_to_grid_dist_matrix {
		float t1 <- machine_time;
		list<string> geoid_list <- blockgroup collect each.geoid;
		list<string> gridid_list <- landuse collect each.id;
		loop i_idx from:0 to:(length(geoid_list)-1) {
			string from_geoid <- geoid_list[i_idx];
			blockgroup from_blockgroup <- first(blockgroup where (each.geoid = from_geoid));
			loop j_idx from:0 to:(length(gridid_list)-1) {
				string to_gridid <- gridid_list[j_idx];
				landuse to_landuse<- first(landuse where (each.id = to_gridid));
				float this_dist <- distance_to(centroid(from_blockgroup), centroid(to_landuse));
				if from_geoid in blockgroup_to_grid_distance_matrix_map.keys {
					blockgroup_to_grid_distance_matrix_map[from_geoid] << to_gridid::this_dist;
				} else {
					blockgroup_to_grid_distance_matrix_map[from_geoid] <- map(to_gridid::this_dist);
				}
			}
		}
		file f <- json_file(blockgroup_to_grid_distance_matrix_json_path, blockgroup_to_grid_distance_matrix_map);
		save f;
		float t2 <- machine_time;
		write "Create block-grid distance matrix: " + string((t2-t1)/1000) + " seconds elasped";
	}
	
	action create_grid_dist_matrix {
		float t1 <- machine_time;
		list<string> gridid_list <- landuse collect each.id;
		loop i_idx from:0 to:(length(gridid_list)-2) {
			string from_gridid <- gridid_list[i_idx];
			landuse from_landuse <- first(landuse where (each.id = from_gridid));
			if from_gridid in grid_distance_matrix_map.keys {
				grid_distance_matrix_map[from_gridid] << from_gridid::0.0;
			} else {
				grid_distance_matrix_map[from_gridid] <- map(from_gridid::0.0);
			}
			loop j_idx from:(i_idx+1) to:(length(gridid_list)-1) {
				string to_gridid <- gridid_list[j_idx];
				landuse to_landuse <- first(landuse where (each.id = to_gridid));
				float this_dist <- distance_to(centroid(from_landuse), centroid(to_landuse));
				grid_distance_matrix_map[from_gridid] << to_gridid::this_dist;
				if to_gridid in grid_distance_matrix_map.keys {
					grid_distance_matrix_map[to_gridid] << from_gridid::this_dist;
				} else {
					grid_distance_matrix_map[to_gridid] <- map(from_gridid::this_dist);
				}
			}
		}
		string last_gridid <- gridid_list[length(gridid_list)-1];
		grid_distance_matrix_map[last_gridid] << last_gridid::0;
		file f <- json_file(grid_distance_matrix_json_path, grid_distance_matrix_map);
		save f;
		float t2 <- machine_time;
		write "Create grid-grid distance matrix: " + string((t2-t1)/1000) + " seconds elasped";
	}

}


species government {
	float finance <- 10.0;
	float revene;
	float expenditure;
	init {
	}
}

species developer {
	// only consider the balance between expenditure for new constructions and subsidies and income of "additional rent".
	// additional rent is the difference between current total rent (calculated by rent_base, i.e. do not consider subsidy here) and total rent at start.
	// this is based on the consideration that total rent at start would be entirely used as fixed cost (maintainment, salary, etc.) and fixed expected profit,
	// and thus additional rent is the superprofits to stimulate the developer to make some changes. 
	float finance <- 0.0;
	float revene_total <- 0.0;
	float expenditure_total <- 0.0;
	float subsidy_total <- 0.0;
	float total_rent_at_start <- 0.0;
	map<int, map<string, float>> revene_detail <- map([]);
	map<int, map<string, float>> expenditure_detail <- map([]);
	init {
		loop c from:0 to:max_cycle+1 {
			revene_detail << c::map(['rent'::0]);
			expenditure_detail << c::map(['construction'::0, 'subsidy'::0]);
		}
	}
	action collect_additional_rent_and_pay_subsidy {
//		float rent_in_this_cycle <- sum((people where each.live_in_kendall) collect (each.represented_nb_people * each.myrent));
		float total_rent_base_this_cycle <- sum((landuse where (each.usage='R')) collect (each.rent_base * each.crt_total_pop));
		float addtional_rent_this_cycle <- (total_rent_base_this_cycle - total_rent_at_start);
		float subsidy_this_cycle <- sum((people where each.live_in_kendall) collect ((each.home_grid.rent_base - each.myrent) * each.represented_nb_people));
		subsidy_this_cycle <- subsidy_this_cycle * 3;
		addtional_rent_this_cycle <- addtional_rent_this_cycle * 3;
		if incentive_policy {
//			subsidy_this_cycle <- addtional_rent_this_cycle * 0.9 * (1 - rent_discount_ratio_low_inc*rent_discount_ratio_less_commuting*rent_discount_ratio_small_scale);
			subsidy_this_cycle <- subsidy_this_cycle * 1.75;
		}
//		write "subsidies = " + subsidy_this_cycle;
//		write "additional rents = " + addtional_rent_this_cycle;
		subsidy_total <- subsidy_total + subsidy_this_cycle;
		expenditure_total <- expenditure_total + subsidy_this_cycle;
		revene_detail[cycle]['rent'] <- addtional_rent_this_cycle;
		revene_total <- revene_total + addtional_rent_this_cycle;
		finance <- finance + addtional_rent_this_cycle - subsidy_this_cycle;
	}
	
	action collect_additional_rent_and_pay_subsidy_new {
		list<people> people_to_collect_rent;
//		float rent_this_cycle <- sum(people_to_collect_rent collect (each.represented_nb_people) * each.home_grid.rent_base);
	}
	
	action potential_eval {
		
	}
}


species blockgroup {
	string geoid;
	float myarea;
	int init_low_inc_pop;
	int init_high_inc_pop;
	float small_size_apt_ratio;
	float rent;
	map<string, float> rent_subgroups;  //only apply for kendall virtual block
	int crt_nb_available_bedrooms;
	int crt_low_inc_pop;
	int crt_high_inc_pop;
	int crt_total_pop;
	float crt_mean_income;
	float crt_pop_density;
	float crt_diversity;
	list<people> attached_residents <- [];
	list<people> attached_workers <- [];
	map<string,int> nb_workers_in_kendall <- ['low'::0, 'high'::0];
	float pop_density;
	bool is_abstract_whole_kendall <- false;
	
	
	action update_current_population (bool quick <- false) {
		if quick = false {
			crt_low_inc_pop <- sum((attached_residents where (each.income=0)) collect each.represented_nb_people);
			crt_high_inc_pop <- sum((attached_residents where (each.income=1)) collect each.represented_nb_people);
			crt_total_pop <- crt_low_inc_pop + crt_high_inc_pop;
			if crt_total_pop > 0 {
				crt_pop_density <- crt_total_pop / (myarea/10000);
				crt_mean_income <- (crt_low_inc_pop*0 + crt_high_inc_pop*1) / crt_total_pop;
			} else {
				crt_pop_density <- 0.0;
				crt_mean_income <- 0.0;
			}
		} else {
			crt_low_inc_pop <- sum(attached_residents collect each.represented_nb_people);
		}
		
		if self = kendall_virtual_block { // update kendall virtual block
			crt_nb_available_bedrooms <- sum(landuse collect each.crt_nb_available_bedrooms);
			do get_rent_subgroups_for_kendall;
			if quick = false {
				if init_low_inc_pop + init_high_inc_pop > 0 {
					small_size_apt_ratio <- init_low_inc_pop / (init_low_inc_pop + init_high_inc_pop);
				} else {
					small_size_apt_ratio <- 0.0;
				}
				if crt_total_pop > 0{
					crt_diversity <- -(crt_low_inc_pop/crt_total_pop) * ln(crt_low_inc_pop/crt_total_pop) - (crt_high_inc_pop/crt_total_pop) * ln(crt_high_inc_pop/crt_total_pop);
				} else {
					crt_diversity <- 0.0;
				}
			}
		}
	}
	
	
	list<people> create_people_given_income_and_blockgroup (int nb_people, int this_income, 
		bool this_in_kendall, int this_represented_nb_people, rgb this_color
	) {
		/*
		 * 
		 * Arguments
		 * ------------------------------------------
		 * nb_people: number of people to be created
		 * this_income: incom type of these people: 0 = low income, 1 = high income
		 * this_in_kendall: whether or not these people live in a blockgroup whose geoid is in the list of Kendall geoids
		 * this_represented_nb_people: how many people are represented by a single agent
		 * this_color: used in aspect, the color of the agent circle
		 */
		 float this_show_size <- 5 + (20-5)*(this_represented_nb_people-1)/(50-1);
		 create people number:nb_people returns:this_population{
				income <- this_income;
				if this_in_kendall {
					list<landuse> home_grid_candidates;
					if this_income = 0 { //give priority to affordable residences for low income people
						home_grid_candidates <- (landuse overlapping myself) where (each.is_affordable and each.crt_nb_available_bedrooms >= this_represented_nb_people);
					} else {  //give priority to normal residences for high income people
						home_grid_candidates <- (landuse overlapping myself) where (each.is_affordable=false and each.crt_nb_available_bedrooms >= this_represented_nb_people);
					}
					if length(home_grid_candidates) = 0 { //accept the other kind of residence if the preferred one is unavailable
						home_grid_candidates <- (landuse overlapping myself) where (each.crt_nb_available_bedrooms >= this_represented_nb_people);
					}
					if length(home_grid_candidates) > 0 {
						home_grid <- one_of(home_grid_candidates);
					} else {
						home_grid <- (landuse where (each.crt_nb_available_bedrooms >= this_represented_nb_people)) closest_to myself;
//						write "Init people warning: a person in kendall block (geoid=" + myself.geoid + ") can not find a grid inside his block with vacant residence, use close available grid instead";
					}
					if home_grid != nil {
						location <- any_location_in(home_grid);
						home_grid.attached_residents << self;
						home_grid.crt_total_pop <- home_grid.crt_total_pop + this_represented_nb_people;
						if this_income = 0 {
							home_grid.crt_low_inc_pop <- home_grid.crt_low_inc_pop + this_represented_nb_people;
							my_rent_type <- 'low_income';
						} else {
							home_grid.crt_high_inc_pop <- home_grid.crt_high_inc_pop + this_represented_nb_people;
							my_rent_type <- 'all_others';
						}
						home_grid.crt_nb_available_bedrooms <- home_grid.crt_nb_available_bedrooms - this_represented_nb_people;
						myrent <- home_grid.rent_subgroups[my_rent_type];
					} else {
						location <- any_location_in(myself);
						write "Init people warning: all landuse grids are fully occupied, the person has to find a location in nonKendall blocks";
					}
						
				} else {
					location <- any_location_in(myself);
					myrent <- myself.rent;
				}
				home_loc <- location;
				home_block <- myself;
				myself.attached_residents << self;
				represented_nb_people <- this_represented_nb_people;
				mycolor <- this_color;
				show_size <- this_show_size;
				live_in_kendall <- this_in_kendall;
			}
		return this_population;
	}
	
	action get_rent_subgroups_for_kendall (string method <- 'min') {  // only apply to kendall virtual block
		if crt_nb_available_bedrooms > 0 {
			loop person_type over: ['all_others', 'low_income', 'less_commuting', 'low_income_and_less_commuting'] {
				if method = 'mean' {
					rent_subgroups[person_type] <- sum(landuse collect (each.rent_subgroups[person_type] * each.crt_nb_available_bedrooms)) / 
						crt_nb_available_bedrooms;
				} else if method = 'min' {
					rent_subgroups[person_type] <- min((landuse where (each.crt_nb_available_bedrooms>0)) collect each.rent_subgroups[person_type]);
				}
			}
		} else {
			rent_subgroups <- ['all_others'::0.0, 'low_income'::0.0, 'less_commuting'::0.0, 'low_income_and_less_commuting'::0.0];
		}
	}
	
	aspect base {
		if show_blockgroup and is_abstract_whole_kendall=false{
			draw shape color: rgb(#lightgrey, 40) border:rgb(#black, 40);
		}
	}
}

species people schedules: shuffle(people){
	int represented_nb_people;
	int settled_time;
	float show_size;
	rgb mycolor <- #grey;
	int income;                   //0:low, 1:high
	bool live_in_kendall;
	bool work_in_kendall;
	float commute_distance;
	float commute_distance_init;
	landuse home_grid;
	point home_loc;
	blockgroup home_block;
	landuse work_grid;
	point work_loc;
	blockgroup work_block;
	map<string,float> distance_to_blockgroups_from_my_work_map;
	map<string,float> distance_to_grids_from_my_work_map;
	list<string> nonkendall_geoids_in_3km_from_work;
	float b_move;
	float b_commute_distance;
	float b_large_size;
	float b_price;
	float b_pop_density;
	float b_inc_disparity;
	matrix<float> b_mat;
	float utility <- 999.0;  //only apply to those who live in kendall, here 999 will be viewed as nan;
	float myrent;
	string my_rent_type <- '';   //only apply to those live in kendall: all_others / low_income / less_commuting / low_income_and_less_commuting
	bool less_commuting_if_move_to_kendall <-false;
	init {
		if income=0 {
			b_move <- b_move_low_inc;
			b_commute_distance <- b_commute_distance_low_inc;
			b_large_size <- b_large_size_low_inc;
			b_price <- b_price_low_inc;
			b_pop_density <- b_pop_density_low_inc;
			b_inc_disparity <- b_inc_disparity_low_inc;
		} else if income=1 {
			b_move <- b_move_high_inc;
			b_commute_distance <- b_commute_distance_high_inc;
			b_large_size <- b_large_size_high_inc;
			b_price <- b_price_high_inc;
			b_pop_density <- b_pop_density_high_inc;
			b_inc_disparity <- b_inc_disparity_high_inc;
		}
		b_mat <- reverse(matrix([b_move, b_commute_distance, b_large_size, b_price, b_pop_density, b_inc_disparity]));
	}
	
	action select_residence (bool must_move <- false){
		float t1 <- machine_time;
		// adjust preference parameters for kendall residents
		if income = 0 {b_move <- live_in_kendall ? b_move_low_inc_kendall : b_move_low_inc;} 
		else {b_move <- live_in_kendall ? b_move_high_inc_kendall : b_move_high_inc;}
		b_mat[0] <- b_move;
		list<string> alternative_geoid_list <- [home_block.geoid]; //current blockgroup will always be in the choice set
		bool kendall_in_choice_set <- false;
		if kendall_virtual_block.crt_nb_available_bedrooms > 0 {  // kendall will be in the choice set as long as there are vacant houses 
			if (distance_to_blockgroups_from_my_work_map['250173523001'] < 3000 and flip(0.2)) or 
			   (distance_to_blockgroups_from_my_work_map['250173523001'] >= 3000 and flip(2*nb_nonkendall_blocks_in_hom_loc_choice/length(nonkendall_geoid_list+1))
			) {
				kendall_in_choice_set <- true;
				alternative_geoid_list << 'kendall';
			}
		}
		// add non-kendall blockgroups as alternatives:
		// 1. randomly select up-to-5 blockgroups from non-kendall blockgroups within 3km from work and with at least 1 vacatn bedroom.
		// 2. then, randomly select up-to-5 blockgroups from all other non-kendall blockgroups with at least 1 vacant bedroom
		list<string> candidate_nearby_nonkendall_blockgroups_with_vacant_residences;
		list<string> candidate_other_nonkendall_blockgroups_with_vacant_residences;
		if consider_nonkendall_resident_capacity {
			 candidate_nearby_nonkendall_blockgroups_with_vacant_residences <- nonkendall_geoids_in_3km_from_work where (
				blockgroup_lookup_map[each].crt_nb_available_bedrooms > 0);
		} else {
			candidate_nearby_nonkendall_blockgroups_with_vacant_residences <- nonkendall_geoids_in_3km_from_work where (
				blockgroup_lookup_map[each].init_low_inc_pop + blockgroup_lookup_map[each].init_high_inc_pop>0
			);
		}
		alternative_geoid_list <-  alternative_geoid_list + 
			(nb_nonkendall_blocks_in_hom_loc_choice among candidate_nearby_nonkendall_blockgroups_with_vacant_residences); 
		
		if consider_nonkendall_resident_capacity {
			candidate_other_nonkendall_blockgroups_with_vacant_residences<- (nonkendall_geoid_list - alternative_geoid_list) where (
				blockgroup_lookup_map[each].crt_nb_available_bedrooms > 0);
		} else {
			candidate_other_nonkendall_blockgroups_with_vacant_residences <- (nonkendall_geoid_list - alternative_geoid_list) where (
				blockgroup_lookup_map[each].init_low_inc_pop + blockgroup_lookup_map[each].init_high_inc_pop>0
			);
		}
		 
		alternative_geoid_list <- alternative_geoid_list + 
			(nb_nonkendall_blocks_in_hom_loc_choice among candidate_other_nonkendall_blockgroups_with_vacant_residences);
		
		// setup utility and make choice
		list<float> move_or_stay_to_these_blocks <- alternative_geoid_list collect 1.0;    // 0=stay, 1=move
		move_or_stay_to_these_blocks[0] <- 0.0;
		list<float> commute_dist_to_these_blocks <- alternative_geoid_list collect (distance_to_blockgroups_from_my_work_map[(each='kendall' ? '250173523001' : each)] / 1000);
		list<float> size_in_these_blocks <- alternative_geoid_list collect (float(flip(blockgroup_lookup_map[each].small_size_apt_ratio)));
		list<float> rent_in_these_blocks <- alternative_geoid_list collect (blockgroup_lookup_map[each].rent);
		if rent_in_these_blocks[0] > 10000 {
			rent_in_these_blocks[0] <- (income=0 ? mean_monthly_rent_low_income : mean_monthly_rent_high_income);
		}
//		if kendall_in_choice_set {
//			string rent_type_in_kendall_block <- get_rent_type([commute_dist_to_these_blocks[1]])[0];
//			rent_in_these_blocks[1] <- kendall_virtual_block.rent_subgroups[rent_type_in_kendall_block];
//		}
		list<float> pop_density_in_these_blocks <- alternative_geoid_list collect (blockgroup_lookup_map[each].crt_pop_density);
		list<float> inc_disparity_to_these_blocks <- alternative_geoid_list collect (abs(blockgroup_lookup_map[each].crt_mean_income - income));
		
		matrix<float> tmp_x <- matrix([move_or_stay_to_these_blocks, commute_dist_to_these_blocks, size_in_these_blocks,
								rent_in_these_blocks, pop_density_in_these_blocks, inc_disparity_to_these_blocks
		]);
		list<float> v <- list(tmp_x.b_mat);
		
		//if kendall in choice set, will use the max v of all grids for him/her to represent kendall
		list<landuse> alternative_grid_list <- [];
		list<float> v_grids <- [];
		list<string> rent_type_in_these_grids;
		if kendall_in_choice_set {
			alternative_grid_list <- landuse where (each.crt_nb_available_bedrooms > 0);
			list<float> move_or_stay_to_these_grids <- alternative_grid_list collect 1.0;
			list<float> commute_dist_to_these_grids <- alternative_grid_list collect (distance_to_grids_from_my_work_map[each.id] / 1000);
			list<float> size_in_these_grids <- alternative_grid_list collect (float(each.is_affordable));
			rent_type_in_these_grids <- get_rent_type(commute_dist_to_these_grids);
			list<float> rent_in_these_grids <- [];
			loop grid_idx from:0 to:length(alternative_grid_list)-1 {
				landuse this_grid <- alternative_grid_list[grid_idx];
				string rent_type_in_this_grid <- rent_type_in_these_grids[grid_idx];
				rent_in_these_grids << this_grid.rent_subgroups[rent_type_in_this_grid];
			}
			list<float> pop_density_in_these_grids <- alternative_grid_list collect (each.crt_pop_density);
			list<float> inc_disparity_to_these_grids <- alternative_grid_list collect (abs(each.crt_mean_income - income));
			matrix<float> tmp_x_grids <- matrix([move_or_stay_to_these_grids, commute_dist_to_these_grids, size_in_these_grids,
									rent_in_these_grids, pop_density_in_these_grids, inc_disparity_to_these_grids
			]);
			v_grids <- list(tmp_x_grids.b_mat);
			float max_v_grid <- max(v_grids);
			v[1] <- max_v_grid;
			// for debug:
			int idx_argmax_v_grid <- v_grids index_of max_v_grid;
			rent_in_these_blocks[1] <- rent_in_these_grids[idx_argmax_v_grid];
			commute_dist_to_these_blocks[1] <- commute_dist_to_these_grids[idx_argmax_v_grid];
		}
		
		list<float> p;
		if represented_nb_people > lower_nb_represented_persons {// if represented_nb_people is big enough, use the real aggreated prob 
			 p <- v_to_p(v);
			 if must_move { //set p[0](stay)=0 and then rescale p vector
			 	p[0] <- 0;
			 	float sum_p <- sum(p);
			 	if sum_p != 0 {p <- p collect (each/sum_p);} else {p <- p collect (1/length(p));}
			 }
		} else { //if too few people make choice, aggreated prob may cause round error, like round([0.3,0.3,0.4] * 1 person) = [0,0,0]
//			int chosen_idx <- to_choose(v); // to avoid this round error, use simulated prob by Monte Carlo method
			int chosen_idx <- (v index_of max(v)); //to avoid this round error, always choose the alt with max_v
			if must_move and chosen_idx = 0 {
				v[0] <- min(v)-1;
				chosen_idx <- (v index_of max(v));
			}
			p <- v collect 0.0;
			p[chosen_idx] <- 1.0;
		}
		// debug
//		if live_in_kendall {
//			write alternative_geoid_list;
//			write v;
//			write p;
//			write '-----------------\n';
//		}

//		if kendall_in_choice_set {
//			write alternative_geoid_list;
//			write "rent = " + rent_in_these_blocks;
//			write "distance = " + commute_dist_to_these_blocks;
//			write "p = "+p;
//			write "Currnet cd = " + commute_distance + ', cd if move to kendall = ' + distance_to_blockgroups_from_my_work_map['250173523001'];
//			write "\n";
////			write round(represented_nb_people * p[1]);
//		}
			
		loop idx from:0 to:length(p)-1 {
			int nb_people_to_make_this_choice <- round(represented_nb_people * p[idx]);
			if idx = 0 {
				// choose to stay, pass
			} else if kendall_in_choice_set and idx = 1 and nb_people_to_make_this_choice > 0 { 
				// choose to move to kendall, do a secenod choice on sepecific landuse grid
				// use a easy way 
				int tmp_nb_people_to_make_this_choice <- nb_people_to_make_this_choice;
				int grid_idx <- 0;
				list<int> sorted_grid_idx <- (v_grids sort_by -each) collect (v_grids index_of each);
				list<landuse> sorted_grids <- sorted_grid_idx collect alternative_grid_list[each];
				loop while: tmp_nb_people_to_make_this_choice>0 and grid_idx<length(alternative_grid_list) {
					landuse destination_grid <- alternative_grid_list[grid_idx];
					int nb_people_to_move_to_this_grid <- min(destination_grid.crt_nb_available_bedrooms, tmp_nb_people_to_make_this_choice);
					tmp_nb_people_to_make_this_choice <- tmp_nb_people_to_make_this_choice - nb_people_to_move_to_this_grid;
					represented_nb_people <- represented_nb_people - nb_people_to_move_to_this_grid;
					home_block.crt_nb_available_bedrooms <- home_block.crt_nb_available_bedrooms + nb_people_to_move_to_this_grid;
					if home_grid != nil {
						home_grid.crt_nb_available_bedrooms <- home_grid.crt_nb_available_bedrooms + nb_people_to_move_to_this_grid;
						kendall_virtual_block.crt_nb_available_bedrooms <- kendall_virtual_block.crt_nb_available_bedrooms + nb_people_to_move_to_this_grid;
					}
					destination_grid.associated_blockgroup.crt_nb_available_bedrooms <- max( 
						destination_grid.associated_blockgroup.crt_nb_available_bedrooms - nb_people_to_move_to_this_grid, 0
					);
					kendall_virtual_block.crt_nb_available_bedrooms <- max(
						kendall_virtual_block.crt_nb_available_bedrooms - nb_people_to_move_to_this_grid, 0
					);
					destination_grid.crt_nb_available_bedrooms <- max(
						destination_grid.crt_nb_available_bedrooms - nb_people_to_move_to_this_grid, 0
					);
					if live_in_kendall = false {  // moving from other places to kendall
						move_to_kendall_pop[income] <- move_to_kendall_pop[income] + nb_people_to_move_to_this_grid;
					}
					
					list<people> people_list_in_destination_grid_with_same_attribute <- destination_grid.attached_residents where (
						each.income=income and each.work_block=work_block and each.work_grid=work_grid and each.my_rent_type=rent_type_in_these_grids[grid_idx]
					);
					if length(people_list_in_destination_grid_with_same_attribute) > 0 {
						ask first(people_list_in_destination_grid_with_same_attribute) {
							represented_nb_people <- represented_nb_people + nb_people_to_move_to_this_grid;
						}
						move_stat << map(['from'::home_loc, 'to'::first(people_list_in_destination_grid_with_same_attribute).home_loc, 
							'pop'::nb_people_to_move_to_this_grid, 'work_in_kendall'::work_in_kendall, 'live_in_kendall'::true]);
					} else {
						people new_generated_people <- copy_myself(
							destination_grid.associated_blockgroup, destination_grid, nil, nb_people_to_move_to_this_grid, 0
						);
						destination_grid.associated_blockgroup.attached_residents << new_generated_people;
						kendall_virtual_block.attached_residents << new_generated_people;
						destination_grid.attached_residents << new_generated_people;
						move_stat << map(['from'::home_loc, 'to'::new_generated_people.home_loc , 'pop'::nb_people_to_move_to_this_grid,
							'work_in_kendall'::work_in_kendall, 'live_in_kendall'::true]);
					}
					
					grid_idx <- grid_idx + 1;
				}
				
				
//				list<float> p_grids;
//				if nb_people_to_make_this_choice > lower_nb_represented_persons {// if nb_people_to_make_this_choice is big enough, use the real aggreated prob 
//					 p_grids <- v_to_p(v_grids);
//				} else { //if too few people make choice, aggreated prob may cause round error, like round([0.3,0.3,0.4] * 1 person) = [0,0,0]
////					int chosen_idx <- to_choose(v_grids);   // to avoid this round error, use simulated prob by Monte Carlo method
//					int chosen_idx <- v_grids index_of max(v_grids);  // to avoid this roud error, always choose the alt with max v
//					p_grids <- v_grids collect 0.0;
//					p_grids[chosen_idx] <- 1.0;
//				}
//				
//				loop grid_idx from:0 to:length(alternative_grid_list)-1 {
//					int nb_people_to_move_to_this_grid <- round(nb_people_to_make_this_choice * p_grids[grid_idx]);
//					write nb_people_to_move_to_this_grid;
//					if nb_people_to_move_to_this_grid > 0 {
//						landuse destination_grid <- alternative_grid_list[grid_idx];
//						// in case that the actually number of vacant bedrooms is less than the number of people going to move to this grid
//						nb_people_to_move_to_this_grid <- min(nb_people_to_move_to_this_grid, destination_grid.crt_nb_available_bedrooms);
//						represented_nb_people <- represented_nb_people - nb_people_to_move_to_this_grid;
//						home_block.crt_nb_available_bedrooms <- home_block.crt_nb_available_bedrooms + nb_people_to_move_to_this_grid;
//						if home_grid != nil {
//							home_grid.crt_nb_available_bedrooms <- home_grid.crt_nb_available_bedrooms + nb_people_to_move_to_this_grid;
//							kendall_virtual_block.crt_nb_available_bedrooms <- kendall_virtual_block.crt_nb_available_bedrooms + nb_people_to_move_to_this_grid;
//						}
//						destination_grid.associated_blockgroup.crt_nb_available_bedrooms <- max( 
//							destination_grid.associated_blockgroup.crt_nb_available_bedrooms - nb_people_to_move_to_this_grid, 0
//						);
//						kendall_virtual_block.crt_nb_available_bedrooms <- max(
//							kendall_virtual_block.crt_nb_available_bedrooms - nb_people_to_move_to_this_grid, 0
//						);
//						destination_grid.crt_nb_available_bedrooms <- max(
//							destination_grid.crt_nb_available_bedrooms - nb_people_to_move_to_this_grid, 0
//						);
//						if live_in_kendall = false {  // moving from other places to kendall
//							move_to_kendall_pop[income] <- move_to_kendall_pop[income] + nb_people_to_move_to_this_grid;
//						}
//						
//						list<people> people_list_in_destination_grid_with_same_attribute <- destination_grid.attached_residents where (
//							each.income=income and each.work_block=work_block and each.work_grid=work_grid and each.my_rent_type=rent_type_in_these_grids[grid_idx]
//						);
//						if length(people_list_in_destination_grid_with_same_attribute) > 0 {
//							ask first(people_list_in_destination_grid_with_same_attribute) {
//								represented_nb_people <- represented_nb_people + nb_people_to_move_to_this_grid;
//							}
//							move_stat << map(['from'::home_loc, 'to'::first(people_list_in_destination_grid_with_same_attribute).home_loc, 
//								'pop'::nb_people_to_move_to_this_grid, 'work_in_kendall'::work_in_kendall, 'live_in_kendall'::true]);
//						} else {
//							people new_generated_people <- copy_myself(
//								destination_grid.associated_blockgroup, destination_grid, nil, nb_people_to_move_to_this_grid, 0
//							);
//							destination_grid.associated_blockgroup.attached_residents << new_generated_people;
//							kendall_virtual_block.attached_residents << new_generated_people;
//							destination_grid.attached_residents << new_generated_people;
//							move_stat << map(['from'::home_loc, 'to'::new_generated_people.home_loc , 'pop'::nb_people_to_move_to_this_grid,
//								'work_in_kendall'::work_in_kendall, 'live_in_kendall'::true]);
//						}
//					}
//				} 
				
			} else if (idx > 1 or (kendall_in_choice_set=false and idx=1)) and nb_people_to_make_this_choice > 0 { 
				//choose to move to blockgroups other than kendall
				blockgroup destination_block <- blockgroup_lookup_map[alternative_geoid_list[idx]];
				represented_nb_people <- represented_nb_people - nb_people_to_make_this_choice;
				home_block.crt_nb_available_bedrooms <- home_block.crt_nb_available_bedrooms + nb_people_to_make_this_choice;
				destination_block.crt_nb_available_bedrooms <- max(destination_block.crt_nb_available_bedrooms - nb_people_to_make_this_choice, 0);
				if live_in_kendall {  //moving from kendall to outside		
					move_out_of_kendall_pop[income] <- move_out_of_kendall_pop[income] + nb_people_to_make_this_choice;
					home_grid.crt_nb_available_bedrooms <- home_grid.crt_nb_available_bedrooms + nb_people_to_make_this_choice;
					kendall_virtual_block.crt_nb_available_bedrooms <- kendall_virtual_block.crt_nb_available_bedrooms + nb_people_to_make_this_choice;
					
				}
				list<people> people_list_in_destination_block_with_same_attribute <- destination_block.attached_residents where (
					each.income=income and each.work_block=work_block and each.work_grid=work_grid
				);
				if length(people_list_in_destination_block_with_same_attribute) > 0 {
					ask first(people_list_in_destination_block_with_same_attribute) {
						represented_nb_people <- represented_nb_people + nb_people_to_make_this_choice;
					}
					move_stat << map(['from'::home_loc, 'to'::first(people_list_in_destination_block_with_same_attribute).home_loc, 
						'pop'::nb_people_to_make_this_choice, 'work_in_kendall'::work_in_kendall, 'live_in_kendall'::live_in_kendall]);
				} else {
					people new_generated_people <- copy_myself(destination_block, nil, nil, nb_people_to_make_this_choice, 0);
					destination_block.attached_residents << new_generated_people;
					move_stat << map(['from'::home_loc, 'to'::new_generated_people.home_loc , 'pop'::nb_people_to_make_this_choice,
						'work_in_kendall'::work_in_kendall, 'live_in_kendall'::live_in_kendall]);
				}
				
			}
		}
		
		/* to accelerate the calculation, do not make choice one by one, use aggregated probability instead.
		 loop times: int(represented_nb_people / lower_nb_represented_persons) {
			int this_chosen_idx <- to_choose(v);
			if this_chosen_idx = 0 {
				// choose to stay, pass
			} else if this_chosen_idx = 1 { // choose to move to kendall
				represented_nb_people <- represented_nb_people - lower_nb_represented_persons;
			} else {
				represented_nb_people <- represented_nb_people - lower_nb_represented_persons;
				home_block.crt_nb_available_bedrooms <- home_block.crt_nb_available_bedrooms + lower_nb_represented_persons;
				blockgroup destination_block <- blockgroup_lookup_map[alternative_geoid_list[this_chosen_idx]];
				list<people> people_list_in_destination_block_with_same_attribute <- destination_block.attached_residents where (
					each.income=income and each.work_block=work_block and each.work_grid=work_grid
				);
				if length(people_list_in_destination_block_with_same_attribute) > 0{
					ask first(people_list_in_destination_block_with_same_attribute) {
						represented_nb_people <- represented_nb_people + lower_nb_represented_persons;
						home_block.crt_nb_available_bedrooms <- max(home_block.crt_nb_available_bedrooms - lower_nb_represented_persons, 0);
					}
				} else {
					people new_generated_people <- copy_myself(destination_block, nil, nil, lower_nb_represented_persons);
					x_tmp <- x_tmp + 1;
					write x_tmp;
				}
			}
		}
		 */
	}
	
	action select_workplace (float prob_in_kendall) {
		if flip(prob_in_kendall) {
			list<landuse> work_grid_candidates <- landuse where (each.worker_capacity-each.crt_nb_workers >= represented_nb_people);
			if length(work_grid_candidates) > 0 {
				work_grid <- any(work_grid_candidates);
				work_grid.crt_nb_workers <- work_grid.crt_nb_workers + represented_nb_people;
				work_grid.attached_workers << self;
				work_loc <- any_location_in(work_grid);
				work_block <- work_grid.associated_blockgroup;
			} else {
				write "Set workplace for Kendall workers: cannot find a inside workplace";
				work_block <- any(blockgroup);
				work_loc <- any_location_in(work_block);
			}
		} else {
			list<float> dist_to_nonkendall_blockgruops_from_my_home <- nonkendall_blockgroups collect (blockgroup_distance_matrix_map[home_block.geoid][each.geoid] / 1000);
			list<float> exp_v <- dist_to_nonkendall_blockgruops_from_my_home collect exp(each* beta_work_allocation);
			int choice_idx <- rnd_choice(exp_v);
			work_block <- nonkendall_blockgroups[choice_idx];
			work_loc <- any_location_in(work_block);
		}
		distance_to_blockgroups_from_my_work_map <- blockgroup_distance_matrix_map[work_block.geoid];
		nonkendall_geoids_in_3km_from_work <- nonkendall_geoid_list where (distance_to_blockgroups_from_my_work_map[each] <= 3000.0);
		if work_grid != nil {
			distance_to_grids_from_my_work_map <- grid_distance_matrix_map[work_grid.id];
		} else {
			distance_to_grids_from_my_work_map <- blockgroup_to_grid_distance_matrix_map[work_block.geoid];	
		}
		if work_grid != nil and home_grid != nil {
			commute_distance <- grid_distance_matrix_map[home_grid.id][work_grid.id];}
		else if work_grid = nil and home_grid = nil {commute_distance <- blockgroup_distance_matrix_map[home_block.geoid][work_block.geoid];}
		else if work_grid != nil {commute_distance <- blockgroup_to_grid_distance_matrix_map[home_block.geoid][work_grid.id];}
		else if home_grid != nil {commute_distance <- blockgroup_to_grid_distance_matrix_map[work_block.geoid][home_grid.id];}	
	}
	
	action calculate_utility {
		if live_in_kendall {
			utility <- commute_distance*b_commute_distance/1000 + b_large_size*(float(home_grid.scale!='S')) + 
				b_price*myrent + b_pop_density*home_grid.crt_pop_density + b_inc_disparity*abs(home_grid.crt_mean_income-income);
		}
	}
	
	list<float> v_to_ev (list<float> v) {
		// calc exp(v) given v, and prevent over flow
		float mean_v <- mean(v);
		list<float> v2 <- v collect (each-mean_v <= 600 ? (each-mean_v) : 600);
		list<float> ev <- v2 collect (exp(each));
		return ev;
	}
	
	list<float> v_to_p (list<float> v) {
		// calc prob given utility (v)
		list<float> ev <- v_to_ev(v);
		float sum_ev <- sum(ev);
		list<float> p <- ev collect (each / sum_ev);
		return p;
	}
	
	int to_choose (list<float> v) {
		// make choice according to utility (v), and return the chosen index
		list<float> ev <- v_to_ev(v);
		int chosen_idx <- rnd_choice(ev);
//		if chosen_idx != 0 {
//			write v;
//			write ev;
//			write chosen_idx;
//			write '\n\n';
//		}
		return chosen_idx;
	}
	
	list<string> get_rent_type (
		list<float> commute_distance_to_alts, float min_commute_distance_before <- min_commute_distance_before_criteria,
		float min_decrease_ratio <- min_decrease_ratio_criteria, float min_decrease_distance <- min_decrease_distance_criteria
	) {
		/* make judegments on what type of rent this person should pay if choose the corresponding residence, mainly decided by chnage of commute distance
		 * there are 3 conditions (A&B&C), all of them must be meet to define the move as "less commuting"
		 * A: the commute distance before must greater than or equal to min_commute_distance_before, 0 = no requirement
		 * B: the decrease of commute distance must greater than or equal to min_decrease_distance, 0 = no requirement
		 * C: (the decrease of commute distance) / (the commute distance before) must be greater than or equal to min_decrease_ratio
		 */
		list<bool> less_commuting_list;
		float commute_distance_km <- commute_distance / 1000;
		if commute_distance_km = 0 {
			less_commuting_list <- commute_distance_to_alts collect false;
		} else {
			less_commuting_list <- commute_distance_to_alts collect (
				(commute_distance_km >= min_commute_distance_before) and
				((commute_distance_km - each) >=min_decrease_distance) and 
				((commute_distance_km - each) / commute_distance_km >= min_decrease_ratio)
			);
		}
		list<string> rent_type <- [];
		loop less_commuting over: less_commuting_list {
			if income=0 and less_commuting {rent_type <- rent_type + 'low_income_and_less_commuting';}
			else if income=1 and less_commuting {rent_type <- rent_type + 'less_commuting';}
			else if income=0 and less_commuting=false {rent_type <- rent_type + 'low_income';}
			else if income=1 and less_commuting=false {rent_type <- rent_type + 'all_others';}
		}
		return rent_type;
	}
	
	people copy_myself (blockgroup new_home_block, landuse new_home_grid, point new_home_loc, int new_represented_nb_people, int new_settled_time) {
		create people number: 1 returns: copied_people {
			represented_nb_people <- (new_represented_nb_people = -1) ? myself.represented_nb_people : new_represented_nb_people;
			show_size <- 5 + (20-5)*(represented_nb_people-1)/(50-1);
			mycolor <- myself.mycolor;
			income <- myself.income;
			live_in_kendall <- (new_home_block in kendall_blockgroups) or (new_home_block = 'kendall');
			if live_in_kendall and new_home_grid = nil {
				// todo: should give this guy a new_home_grid, or set new_home_block out of kendall
				write "Error: new_home_grid should be provided if new_home_block is in kendall";
			}
			work_in_kendall <- myself.work_in_kendall;
			home_block <- new_home_block;
			home_grid <- new_home_grid;
			if new_home_loc = nil {
				if home_grid != nil {
					home_loc <- any_location_in(home_grid);
				} else {
					home_loc <- any_location_in(home_block);
				}
			} else {
				home_loc <- new_home_loc;
			}
			settled_time <- (new_settled_time = -1) ? myself.settled_time : new_settled_time;
			work_grid <- myself.work_grid;
			work_loc <- myself.work_loc;
			work_block <- myself.work_block;
			commute_distance <- myself.commute_distance;   //just for calculation of my_rent_type: whether commute_distance decrease
			commute_distance_init <- myself.commute_distance_init;
			float new_commute_distance;
			if work_grid != nil and home_grid != nil {new_commute_distance <- grid_distance_matrix_map[home_grid.id][work_grid.id];}
			else if work_grid = nil and home_grid = nil {new_commute_distance <- blockgroup_distance_matrix_map[home_block.geoid][work_block.geoid];}
			else if work_grid != nil {new_commute_distance <- blockgroup_to_grid_distance_matrix_map[home_block.geoid][work_grid.id];}
			else if home_grid != nil {new_commute_distance <- blockgroup_to_grid_distance_matrix_map[work_block.geoid][home_grid.id];}
			if home_grid != nil {
				my_rent_type <- get_rent_type([new_commute_distance/1000])[0];
//				// debug	
//				write "after copy: income = " + income + ", kendall rent type: from " + commute_distance/1000 + " km to " + new_commute_distance/1000 + ' km, type=' + my_rent_type;
				myrent <- home_grid.rent_subgroups[my_rent_type];
			} else {
				my_rent_type <- '';
				myrent <- home_block.rent;
			}
			commute_distance <- new_commute_distance;
			distance_to_blockgroups_from_my_work_map <- myself.distance_to_blockgroups_from_my_work_map;
			distance_to_grids_from_my_work_map <- myself.distance_to_grids_from_my_work_map;
			nonkendall_geoids_in_3km_from_work <- myself.nonkendall_geoids_in_3km_from_work;
			b_move <- myself.b_move;
			b_commute_distance <- myself.b_commute_distance;
			b_large_size <- myself.b_large_size;
			b_price <- myself.b_price;
			b_pop_density <- myself.b_pop_density;
			b_inc_disparity <- myself.b_inc_disparity;
			b_mat <- myself.b_mat;
		}
		return first(copied_people);
	}
	
	action update_landuse_and_blockgroup_population_info (bool only_kendall <- false, bool quick <- false) {
		if quick {
			ask landuse where (each.usage='R') {do update_current_population (true);}
			if only_kendall {
				ask kendall_virtual_block {do update_current_population (true);}
			} else {
				ask blockgroup {do update_current_population (true);}
			}
		} else {
			ask landuse {do update_current_population (false);}
			if only_kendall {
				ask kendall_virtual_block {do update_current_population (false);}
			} else {
				ask blockgroup {do update_current_population (false);}
			}
		}
	}
	
	aspect base {
		if show_workplace {
			location <- work_loc;
		} else {
			location <- home_loc;
		}
		if show_people = 'Show all people' or (show_people = 'Show Kendall people only' and 
			((show_workplace and work_grid != nil) or (!show_workplace and home_grid != nil))
		) {
			if show_size > 0 {draw circle(show_size) color: mycolor;}
		} 
		if show_commute {
			float commute_line_width;
			rgb commute_line_color;
			if home_grid != nil and work_grid != nil {
				commute_line_width <- 15.0;
				commute_line_color <- rgb(#black, 120);
			} else if (home_grid != nil or work_grid != nil) {
				if is_focused {
					commute_line_width <- 1.0;
					commute_line_color <- rgb(#grey, 70);
				} else {
					commute_line_width <- 2.0;
					commute_line_color <- rgb(#grey, 70);
				}
				
			} else {
				if is_focused {
					commute_line_width <- 1.0;
					commute_line_color <- rgb(#grey, 20);
				} else {
					commute_line_width <- 2.0;
					commute_line_color <- rgb(#grey, 70);
				}
			}
			if (show_people = 'Show all people') or (show_people = 'Show Kendall people only' and show_workplace and work_grid != nil) 
				or (show_people = 'Show Kendall people only' and show_workplace = false and home_grid != nil) {
					if show_size > 0 {
						draw polyline([home_loc, work_loc]) color:commute_line_color width:commute_line_width end_arrow:25;
					}
			}
		} 

	}
}

species landuse {
	string id;
	string usage;
	float far;
	int max_height;
	string scale;
	float myarea <- 6400.0;
	float base_area;
	float rent_base;
	map<string, float> rent_discount_ratio <- map(['all'::1.0, 'low_income'::1.0, 'small_scale'::1.0, 'less_commuting'::1.0]);
	map<string, float> rent_shift <- map(['all'::0.0, 'low_income'::0.0, 'small_scale'::0.0, 'less_commuting'::0.0]);
	map<string, float> rent_subgroups;
//	float rent update:max(rent_base*rent_discount_ratio + rent_shift, 0.0);
	int init_vacant_nb_bedrooms;  // read from appartment shapefile, maynot compatible with population and residence-landuse information?
	rgb mycolor;
	int resident_capacity <- 0;
	int crt_high_inc_pop <- 0;
	int crt_low_inc_pop <- 0;
	int crt_total_pop <- 0;
	int crt_nb_available_bedrooms;
	float crt_pop_density;
	float crt_mean_income;
	list<people> attached_residents <- [];
	int worker_capacity <- 0;
	int crt_nb_workers <- 0;
	list<people> attached_workers;
	bool is_affordable;
	blockgroup associated_blockgroup;
	list<landuse> my_grid_neighors;
	map<string,float> crt_local_finance;
	map<string,int> target_subgroups <- ['low_income'::0, 'less_commuting'::0];
	float crt_potential <- 0.0;
	list<blockgroup> object_blocks;
	init {
		if usage = 'R' {
			is_affordable <- scale='S';
			resident_capacity <- int(myarea * far / (is_affordable ? residence_area_small_scale : residence_area_large_scale));
			mycolor <- is_affordable ? #red : #orange;
		}
		if usage = 'O' {
			worker_capacity <- int(myarea * far / (scale='S' ? office_area_small_scale : office_area_large_scale));
			mycolor <- #blue;
		}
		if usage = 'vacant' {
			mycolor <- #grey;
		}
		my_grid_neighors <- landuse where (each.id in grid_neighbors_map[id]);
		crt_nb_available_bedrooms <- resident_capacity;
		object_blocks <- nonkendall_blockgroups where (
			blockgroup_to_grid_distance_matrix_map[each.geoid][self.id] <= distance_from_object_block_to_grid_km*1000
		);
	}
	
	action get_object_blockgroups (float distance_km) {
		object_blocks <- nonkendall_blockgroups where (
			blockgroup_to_grid_distance_matrix_map[each.geoid][self.id] <= distance_km*1000
		);
//		write "id="+id + ', length(object_blocks)='+length(object_blocks);
	}
	
	action get_associated_blockgroup {
		associated_blockgroup <- first(blockgroup where (each covers self));
		if associated_blockgroup = nil {
			associated_blockgroup <- first(blockgroup where (each covers centroid(self)));
		}
	}
	
	action update_current_population (bool quick <- false) {
		if quick = false {
			crt_low_inc_pop <- sum((attached_residents where (each.income=0)) collect each.represented_nb_people);
			crt_high_inc_pop <- sum((attached_residents where (each.income=1)) collect each.represented_nb_people);
			crt_total_pop <- crt_low_inc_pop + crt_high_inc_pop;
		} else {
			crt_total_pop <- sum(attached_residents collect each.represented_nb_people);
		}
		crt_nb_available_bedrooms <- max(resident_capacity - crt_total_pop, 0);
		if quick = false {
			int total_pop_around <- sum(my_grid_neighors collect each.crt_total_pop);
			int total_high_inc_pop_around <- sum(my_grid_neighors collect each.crt_high_inc_pop);
			crt_pop_density <- total_pop_around / (0.64 * length(my_grid_neighors));
			if total_pop_around > 0 {
				crt_mean_income <- total_high_inc_pop_around / total_pop_around;
			} else {
				crt_mean_income <- 0.0;
			}
		}
	}
	
	float new_constructions (float far_multiplier, float far_shift, bool add_to_developer_finance <- true) {
		float far_before <- far;
		float construction_fee_multiplier <- 6.0;
		far <- max(0, far*far_multiplier+far_shift);
		float add_floor_area <- 0.0;
		float construction_cost <- 0.0;
		if far > far_before {
			add_floor_area <- (far - far_before) * myarea;
		}
		if add_floor_area > 0 and add_floor_area < 500 {
			construction_cost <- 270 * add_floor_area * construction_fee_multiplier;
		} else if add_floor_area >= 500 and add_floor_area < 2000 {
			construction_cost <- 240 * add_floor_area * construction_fee_multiplier;
		} else if add_floor_area >= 2000 and add_floor_area < 4000 {
			construction_cost <- 200 * add_floor_area * construction_fee_multiplier;
		} else if add_floor_area >= 4000 {
			construction_cost <- 160 * add_floor_area * construction_fee_multiplier;
		}
		if add_to_developer_finance {
			ask developer {
				expenditure_total <- expenditure_total + construction_cost;
				finance <- finance - construction_cost;
				expenditure_detail[cycle]['construction'] <- expenditure_detail[cycle]['construction'] + construction_cost;
			}
		}
		return construction_cost;
	}
	
	action change_landuse (string new_usage, string new_scale <- 'same'){
		string usage_before <- usage;
		string scale_before <- scale;
		if new_usage != 'same' {usage <- new_usage;}
		if new_scale != 'same' {scale <- new_scale;}
		worker_capacity <- 0;
		resident_capacity <- 0;	
		list<people> people_who_need_to_find_new_workplace <- [];
		list<people> people_who_need_to_find_new_residence <- [];
		if usage = 'R' {
			is_affordable <- scale='S';
			resident_capacity <- int(myarea * far / (is_affordable ? residence_area_small_scale : residence_area_large_scale));
			mycolor <- is_affordable ? #red : #orange;
			do update_current_population;
			if usage_before = 'O' {
				people_who_need_to_find_new_workplace <- people_who_need_to_find_new_workplace + attached_workers;
			}
			attached_workers <- [];
			crt_nb_workers <- 0;
			if usage_before = 'R' and crt_total_pop > resident_capacity {
				int crt_total_pop_copy <- crt_total_pop;  //can not reduce crt_total_pop directly, as select_residence will do this again, so make a copy;
				loop while: crt_total_pop_copy > resident_capacity {
					people this_move_out_people <- one_of(attached_residents);
					people_who_need_to_find_new_residence << this_move_out_people;
					crt_total_pop_copy <- crt_total_pop_copy - this_move_out_people.represented_nb_people;
					remove this_move_out_people from: attached_residents;
				}
			}
		}
		if usage = 'O' {
			worker_capacity <- int(myarea * far / (scale='S' ? office_area_small_scale : office_area_large_scale));
			mycolor <- #blue;
			if usage_before = 'R' {
				people_who_need_to_find_new_residence <- people_who_need_to_find_new_residence + attached_residents;
			}
			attached_residents <- [];
			if usage_before = 'O' and crt_nb_workers > worker_capacity {
				loop while: crt_nb_workers > worker_capacity {
					people this_find_new_workplace_people <- one_of(attached_workers);
					people_who_need_to_find_new_workplace << this_find_new_workplace_people;
					crt_nb_workers <- crt_nb_workers - this_find_new_workplace_people.represented_nb_people;
					remove this_find_new_workplace_people from: attached_workers;
				}
			}
		}
		if usage = 'vacant' {
			mycolor <- #grey;
		}
		ask people_who_need_to_find_new_workplace {
			do select_workplace(0.25);
		}
		ask people_who_need_to_find_new_residence {
			do select_residence(true);
		}
	
	}
	
	action set_rent_policy (map new_rent_discount_ratio<-["all"::1.0], map new_rent_shift <- ["all"::0.0]) {
		/*
		 * new_rent_discount_ratio / new_rent_shift are both maps, keys="all", "low_income", "small_scale", "less_commuting"
		 * if some keys are not set, use the current value
		 */
		 loop ratio_key over: new_rent_discount_ratio.keys {
		 	if ratio_key in rent_discount_ratio.keys {rent_discount_ratio[ratio_key] <- new_rent_discount_ratio[ratio_key];}
		 }
		 loop shift_key over:  new_rent_shift.keys {
		 	if shift_key in rent_shift.keys {rent_shift[shift_key] <- new_rent_shift[shift_key];}
		 }
		 float small_scale_rent_discount_ratio <- (scale='S') ?  rent_discount_ratio['small_scale'] : 1.0;
		 float small_scale_rent_shift <- (scale='S') ? rent_shift['small_scale'] : 0.0;
		 rent_subgroups['low_income'] <- max(
		 	0,
		 	rent_base * rent_discount_ratio['all'] * small_scale_rent_discount_ratio * rent_discount_ratio['low_income']
		 		      + rent_shift['all'] + small_scale_rent_shift + rent_shift['low_income']
		 );
		 rent_subgroups['less_commuting'] <- max(
		 	0,
		 	rent_base * rent_discount_ratio['all'] * small_scale_rent_discount_ratio * rent_discount_ratio['less_commuting']
		 		      + rent_shift['all'] + small_scale_rent_shift + rent_shift['less_commuting']
		 );
		 rent_subgroups['low_income_and_less_commuting'] <- max(
		 	0,
		 	rent_base * rent_discount_ratio['all'] * small_scale_rent_discount_ratio * rent_discount_ratio['low_income'] * rent_discount_ratio['less_commuting']
		 		      + rent_shift['all'] + small_scale_rent_shift + rent_shift['low_income'] + rent_shift['less_commuting']
		 );
		 rent_subgroups['all_others'] <- max(
		 	0, 
		 	rent_base * rent_discount_ratio['all'] * small_scale_rent_discount_ratio
		 	          + rent_shift['all'] + small_scale_rent_shift
		 );
		
	}
	
	action update_attached_resident_rent {
		ask attached_residents where (each.settled_time<=1) {
			myrent <- myself.rent_subgroups[my_rent_type];
		}
	}
	
	action target_subgroups_calc (float distance_km <- 3.0, map<string,float> rent_type_weights <- ['low_income'::0.0, 'less_commuting'::0.0]) {
		if usage = 'O' {
			target_subgroups <- ['low_income'::0, 'less_commuting'::0];
		} else {
			do get_object_blockgroups(distance_km);
			loop this_rent_type over:rent_type_weights.keys {
				if this_rent_type = 'low_income' and rent_type_weights[this_rent_type] > 0 {
					list<people> all_residents_in_object_blocks <- object_blocks accumulate each.attached_residents;
					int nb_low_inc_residents_in_object_block <- sum(
						(all_residents_in_object_blocks where (each.income=0 and each.live_in_kendall=false and each.settled_time>=4)) 
						collect each.represented_nb_people
					);
					target_subgroups['low_income'] <- nb_low_inc_residents_in_object_block;
//					crt_potential <- crt_potential + nb_low_inc_residents_in_object_block * rent_type_weights['low_income'];
				}
				if this_rent_type = 'less_commuting' and rent_type_weights[this_rent_type] > 0 {
					list<people> qualified_workers_in_object_blocks <- [];
					loop this_block over: object_blocks {
						float new_commute_distance <- blockgroup_to_grid_distance_matrix_map[this_block.geoid][self.id];
						list<people> this_qualified_workers <- this_block.attached_workers where (
							each.live_in_kendall=false and each.commute_distance/1000 >= min_commute_distance_before_criteria and
							(each.commute_distance - new_commute_distance) / (each.commute_distance) >= min_decrease_ratio_criteria and
							(each.commute_distance/1000 - new_commute_distance/1000) >= min_decrease_distance_criteria
						);
//						write this_qualified_workers;
						qualified_workers_in_object_blocks <- qualified_workers_in_object_blocks + this_qualified_workers;
					}
					int nb_less_commuting_workers_in_object_block <- sum(qualified_workers_in_object_blocks collect each.represented_nb_people);
					target_subgroups['less_commuting'] <- nb_less_commuting_workers_in_object_block;
//					crt_potential <- crt_potential + nb_less_commuting_workers_in_object_block * rent_type_weights['less_commuting'];
				}		
			}
//			write "target_subgroups = " + target_subgroups;
		}
		
	}
	
	
	action potential_calc(
		map<string,float> rent_type_weights <- ['low_income'::0.0, 'less_commuting'::0.0],
		int max_target_low_income_pop, int min_target_low_income_pop,
		int max_target_less_commuting_pop, int min_target_less_commuting_pop
	) {
		if usage='O' {
			crt_potential <- -1000.0;
		} else {
			crt_potential <- 0.0;
			loop this_rent_type over:rent_type_weights.keys {
				if this_rent_type = 'low_income' and rent_type_weights[this_rent_type] > 0 {
					crt_potential <- crt_potential + rent_type_weights['low_income'] * (target_subgroups['low_income']-min_target_low_income_pop) / (max_target_low_income_pop-min_target_low_income_pop);
				}
				if this_rent_type = 'less_commuting' and rent_type_weights[this_rent_type] > 0 {
					crt_potential <- crt_potential + rent_type_weights['less_commuting'] * (target_subgroups['less_commuting']-min_target_less_commuting_pop) / (max_target_less_commuting_pop-min_target_less_commuting_pop);
				}
			}
		}
	}
	
	aspect base {
		draw shape color: rgb(mycolor, 50) border:rgb(#black, 100);
	}
	
}

species building {
	aspect base {
		if show_building {
			draw shape color:#orange;
		}
	}
}

experiment gui type: gui {
	float seed <- 10.0;
//	float minimum_cycle_duration <- 1#sec;
	parameter "Incentive Policy" var:incentive_policy category:"Incentive Policy";
	parameter "Dynamic Incentive Policy" var:dynamic_policy category:"Incentive Policy" disables:[
//		new_construction_grids_string, new_construction_far_string, new_construction_scale, rent_discount_grids_string, 
		construction_intensity, rent_discount_ratio_all,
		rent_discount_ratio_low_inc, rent_discount_ratio_less_commuting, rent_discount_ratio_small_scale
	] enables: [
		diversity_target, low_inc_pop_ratio_target, commute_distance_decrease_target, building_energy_target
	];
//	parameter "New Construction Grids" var:new_construction_grids_string category:"Incentive Policy";
//	parameter "New Construction FAR" var:new_construction_far_string category:"Incentive Policy";
//	parameter "New Construction Scale" var:new_construction_scale category:"Incentive Policy" among:['Small', 'Large'];
	parameter "Construction Intensity" var:construction_intensity category:"Incentive Policy" min:0.0 max:1.0;
//	parameter "Rent Discount Grids" var:rent_discount_grids_string category:"Incentive Policy";
//	parameter "Rent Discount for All People" var:rent_discount_ratio_all max:1.0 min:0.1 category:"Incentive Policy";
	parameter "Rent Discount for Low Income People" var:rent_discount_ratio_low_inc max:1.0 min:0.5 category:"Incentive Policy";
	parameter "Rent Discount for Less Commuting People" var:rent_discount_ratio_less_commuting max:1.0 min:0.5 category:"Incentive Policy";
	parameter "Rent Discount for Small Scale Residences" var:rent_discount_ratio_small_scale max:1.0 min:0.5 category:"Incentive Policy";
	parameter "Diversity Target" var:diversity_target max:0.7 min:0.62 category:"Incentive Policy" ;
	parameter "Affordability Target" var:low_inc_pop_ratio_target max:0.65 min:0.35 category:"Incentive Policy";
	parameter "Commuting Distance Target" var:commute_distance_decrease_target max:0.6 min:-0.05 category:"Incentive Policy";
	parameter "Building Energy Target" var:building_energy_target max:60.0 min:50.0 category:"Incentive Policy";
//	parameter "Diversity Weight" var:diversity_weight category:"Incentive Policy" ;
//	parameter "Affordability Weight" var:low_inc_pop_ratio_weight category:"Incentive Policy";
//	parameter "Commuting Energy Weight" var:commute_distance_weight category:"Incentive Policy";
//	parameter "Building Energy Weight" var:building_energy_weight category:"Incentive Policy";

	parameter "Show Blockgroup" var:show_blockgroup category:Viz;
	parameter "Focus on Kendall Square" var:focus_on_kendall category:Viz;
	parameter "Show Legend" var:show_legend category:Viz;
	parameter "Show People" var:show_people among: ['Show all people', 'Show Kendall people only', 'Hide all people'] category:Viz;
	parameter "Show Workplace" var:show_workplace category:Viz;
	parameter "Show Commute" var:show_commute category:Viz;
	parameter "Show Move" var:show_move category:Viz;
	
	parameter "Move (vs. Stay)" var:b_move_low_inc category: "Preference Parameters: Low Income People";
	parameter "Commute Distance (km)" var:b_commute_distance_low_inc category: "Preference Parameters: Low Income People";
	parameter "Size: Large (vs. Small)" var:b_large_size_low_inc category: "Preference Parameters: Low Income People";
	parameter "Monthly Rent ($)" var:b_price_low_inc category: "Preference Parameters: Low Income People";
	parameter "Population Density (persons/ha)" var:b_pop_density_low_inc category: "Preference Parameters: Low Income People";
	parameter "Income Disparity" var:b_inc_disparity_low_inc category: "Preference Parameters: Low Income People";
	
	parameter "Move (vs. Stay) " var:b_move_high_inc category: "Preference Parameters: High Income People";
	parameter "Commute Distance (km) " var:b_commute_distance_high_inc category: "Preference Parameters: High Income People";
	parameter "Size: Large (vs. Small) " var:b_large_size_high_inc category: "Preference Parameters: High Income People";
	parameter "Monthly Rent ($) " var:b_price_high_inc category: "Preference Parameters: High Income People";
	parameter "Population Density (persons/ha) " var:b_pop_density_high_inc category: "Preference Parameters: High Income People";
	parameter "Income Disparity " var:b_inc_disparity_high_inc category: "Preference Parameters: High Income People";
	
	
	output {
		display base type:java2D{
			species blockgroup aspect: base;
			species landuse aspect: base;
			species people aspect: base;
			species building aspect:base;
			graphics 'move' {
				if show_move {
					loop this_move over: move_stat {
						if show_people = 'Show all people' {
							draw polyline([this_move['from'], this_move['to']]) color:#red width:int(this_move['pop']) end_arrow:25;
						} else if (show_people = 'Show Kendall people only' and show_workplace=false and this_move['live_in_kendall']) {
							draw polyline([this_move['from'], this_move['to']]) color:#red width:int(this_move['pop']) end_arrow:25;
						} else if (show_people = 'Show Kendall people only' and show_workplace and this_move['work_in_kendall']) {
							draw polyline([this_move['from'], this_move['to']]) color:#red width:int(this_move['pop']) end_arrow:25;
						}
						
					}
				}
			}
			overlay position: {10, 10} size: {200 #px, 220 #px} transparency: (show_legend? 1 : 0)
				border: (show_legend ? #black : nil) background: #white {
            	if(show_legend){
	            	float y <- 20#px;
	            	map<string,rgb> landuse_color <- ['Affordable Residence'::#red, 'Normal Residence'::#orange, 'Office'::#blue, 'Vacant'::#grey];
            	  	loop landuse_type over: landuse_color.keys {
            	  		draw square(25#px) at: { 20#px, y } color: rgb(landuse_color[landuse_type], 50) border: #white;
            	  		draw landuse_type at: { 40#px, y + 8#px } color: #black font: font("SansSerif", 15, #bold);
            	  		y <- y + 35#px;
            	  	}	 
            	  	y <- y + 10#px;
            	  	map<string,rgb> people_color <- ['Low Income People'::#red, 'High Income People'::#blue];
            	  	loop people_type over: people_color.keys {
            	  		draw circle(8#px) at: { 20#px, y } color: people_color[people_type] border: #white;
            	  		draw people_type at: { 40#px, y + 5#px } color: #black font: font("SansSerif", 15, #bold);
            	  		y <- y + 35#px;
            	  	}             
            	}
            }
		}
		display kendall_performance type:java2D{
			chart "Total Population in Kendall" type: series background: #white position: {0.0,0.0} size: {0.33, 0.33} {
				data "All" value: kendall_virtual_block.crt_total_pop color:#green;
				data "Low Income" value:kendall_virtual_block.crt_low_inc_pop color:#red;
				data "High Income" value:kendall_virtual_block.crt_high_inc_pop color:#blue;
			}
//			chart "Moving Population Related with Kendall" type: series background: #white position: {0.33,0.0} size: {0.67,0.33} y_range:{0,crt_move_y_axis_upper_range}{
			chart "Moving Population Related with Kendall" type: series background: #white position: {0.33,0.0} size: {0.67,0.33} y_range:{0,200}{
				data "To Kendall: Low Income" value: move_to_kendall_pop[0] color: #green;
				data "To Kendall: High Income" value: move_to_kendall_pop[1] color: #blue;
				data "Out of Kendall: Low Income" value:move_out_of_kendall_pop[0] color: #red;
				data "Out of Kendall: High Income" value:move_out_of_kendall_pop[1] color: #orange;
			}
			chart "Occupancy" type:series background: #white position:{0.0,0.33} size:{0.33,0.34}{
				data "Overall" value: kendall_occupancy['overall'] color: #green;
//				data "Small Scale" value: kendall_occupancy['small'] color:#red;
//				data "Large Scale" value: kendall_occupancy['large'] color:#blue;
			}
//			chart "Commute Distance" type:series background: #white position:{0.33,0.33} size:{0.34,0.34} {
//				data "All" value: mean_commute_distance['all'] color:#red;
//				data "Work in Kendall" value: mean_commute_distance['work_in_kendall'] color:#green;
//				data "Live in Kendall" value: mean_commute_distance['live_in_kendall'] color:#blue;
////				data "Work or Live in Kendall" value: mean_commute_distance['work_or_live_in_kendall'] color:#yellow;
//			}

			chart "Commute Distance Decrease" type:series background:#white position:{0.33,0.33} size:{0.33,0.34} y_range:{-0.1,0.6}{
				data "Mean" value:commute_distance_decrease['mean'] color:#red;
			}

			chart "Residence Energy" type:series background: #white position:{0.67,0.33} size:{0.33,0.34} y_range:{50,60}{
				data "Eerngy per Kendall Resident" value: residence_energy_per_person color: #red;
			}
			
//			chart "Population Structure" type:series background: #white position:{0.0,0.67} size:{0.33,0.33}{
//				data "Diversity" value: kendall_diversity color: #red;
//				data "Low Income Proportion" value: kendall_low_inc_ratio color: #blue;
//			}
			
			chart "Affordability" type:series background: #white position:{0.0,0.67} size:{0.33,0.33} y_range:{0.3,0.65}{
				data "Low Income Proportion" value: kendall_low_inc_ratio color: #red;
			}
			
			chart "Diversity" type:series background:#white position:{0.33,0.67} size:{0.34,0.33} y_range:{0.63,0.7} {
				data "Diversity" value: kendall_diversity color: #red;
			}
//			chart "Resident Mean Utility" type:series background: #white position:{0.33,0.67} size:{0.34,0.33} {
//				data "Overall" value: kendall_resident_utility['total']-kendall_resident_utility_at_start['total'] color: #green;
//				data "Low Income" value: kendall_resident_utility['low_inc']-kendall_resident_utility_at_start['low_inc'] color: #red;
//				data "High Income" value: kendall_resident_utility['high_inc']-kendall_resident_utility_at_start['high_inc'] color: #blue;
//			}
			chart "Developer Finance" type:series background: #white position:{0.67,0.67} size:{0.33,0.33}{
				data "Finance" value: the_developer.finance color: #blue;
				data "Expenditure" value: the_developer.expenditure_total color: #red;
				data "Revene" value: the_developer.revene_total color: #green;
			}
		}
//		monitor "Total population" value: sum(people collect (each.represented_nb_people));
//		monitor "Total population moved" value: sum((people where (each.settled_time<1)) collect (each.represented_nb_people));
//		inspect name:'Kendall' value:kendall_virtual_block type:agent;
	}
	
}