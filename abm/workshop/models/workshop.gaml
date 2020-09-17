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
//	file<geometry> buildings_shp <- file<geometry>("../includes/BuildingsLatLongBlock.shp");
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
	float far_discount <- 0.7;    //it seems that current landuse plan provides too many residences
	list<map> move_stat <- [];
	bool consider_nonkendall_resident_capacity <- false;
	
	int lower_nb_represented_persons <- 3;     //mainly used for kendall landuse, and as the unit for a real move
	int upper_nb_represented_persons <- 120;   //mainly used for outside blockgroup
	int nb_nonkendall_blocks_in_hom_loc_choice <- 10;
	float ratio_of_people_considering_move_at_each_step <- 0.05;
	float mean_monthly_rent_low_income <- 1112.0; //25% percentile for all rents, used as rent of current residence for low inc person if block mean rent unavailabe
	float mean_monthly_rent_high_income <- 1429.0; //75% percentile for all rents, used as rent of current residence for high inc person if block mean rent unavailabe
	float beta_work_allocation <- -0.5;
	
	// preference parameters
	float b_move_low_inc <- -1.43;
	float b_commute_distance_low_inc <- -0.88; //-1.28
	float b_large_size_low_inc <- 0.0; //0.24
	float b_price_low_inc <- -0.003;    //-0.004
	float b_pop_density_low_inc <- -0.00; //-0.0056
	float b_inc_disparity_low_inc <- -0.00;  //-0.21
	
	float b_move_high_inc <- -1.29;
	float b_commute_distance_high_inc <- -0.81;  //-1.31
	float b_large_size_high_inc <- 0.0;  //0.52
	float b_price_high_inc <- -0.0006;    //-0.0012
	float b_pop_density_high_inc <- -0.0; //-0.0123
	float b_inc_disparity_high_inc <- -0.0;  //-0.52
	
	//policy parameters
	float construction_cost_per_m2 <- 1000.0;
	int subsidy_low_inc <- 50;
	int sbusidy_less_commute <- 500;
	float commute_distance_goal <- 2.0;
	
	
	// performance indices
	float kendall_occupancy;
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
			rent <- max(rent_base*rent_discount_ratio + rent_shift, 0.0);
		}
		kendall_envelope <- envelope(landuse) * 1;
		
//		create building from:buildings_shp;
		create government number:1;
		create developer number:1;
		
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
			if crt_nb_available_bedrooms > 0 {
				rent <- sum(landuse collect (each.rent * each.crt_nb_available_bedrooms)) / crt_nb_available_bedrooms;
			} else {
				rent <- 0.0;
			}
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
//		write 'negative number of vacant hoess: ' + (landuse where (each.crt_nb_available_bedrooms<0)) collect each.id;
		ask landuse {
			do update_current_population;
		}
		ask blockgroup {
			do update_current_population;
		}
		do kendall_statistics;
	}
	
	reflex main {
		write "\nNew step start:";
		move_to_kendall_pop <- map([0::0, 1::0]);
		move_out_of_kendall_pop <- map([0::0, 1::0]);
		move_stat <- [];
		
		// apply Kenall policy 
		if cycle > 4 {
			write "Landuse rent discount: 0.01";
			ask landuse {
				rent_discount_ratio <- 0.8;
			}
		}
		
		ask people {
			settled_time <- settled_time + 1;   
//			if (work_in_kendall=true or flip(ratio_of_people_considering_move_at_each_step)) and settled_time>=5 {
			if (flip(ratio_of_people_considering_move_at_each_step)) and settled_time>=5 {
				do select_residence;
			}
		}
		write 'people over, vacant house=' + kendall_virtual_block.crt_nb_available_bedrooms;
		ask landuse {
			do update_current_population;
		}
		ask blockgroup {
			do update_current_population;
		}
		write 'update over, vacant house=' + kendall_virtual_block.crt_nb_available_bedrooms;
		do kendall_statistics;
		write 'statistic over, vacant house=' + kendall_virtual_block.crt_nb_available_bedrooms;
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
	
	reflex halting when: cycle>12 {
		do pause;
	}
	
	action reference_info {
//		write machine_time;
//		write 'negative number of vacant hoess: ' + (landuse where (each.crt_nb_available_bedrooms<0)) collect each.id;
		list<people> tmp_people <- people where (each.work_in_kendall and each.live_in_kendall=false and each.represented_nb_people>0);
//		write "Debug info start: \n";
//		write '\n\n==================\nPoeple work in Kendall but live outside: ' + sum(tmp_people collect each.represented_nb_people);
//		write tmp_people;
//		write tmp_people collect each.represented_nb_people;
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
		ask blockgroup where (each.is_abstract_whole_kendall = false) {
			int this_represented_nb_people;
			bool this_in_kendall <- (geoid in kendall_geoid_list) ? true : false;
			list<people> this_low_income_population_work_in_kendall <- [];
			list<people> this_high_income_population_work_in_kendall <- [];
			list<people> this_low_income_population_work_outof_kendall <- [];
			list<people> this_high_income_population_work_outof_kendall <- [];
			
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
				list<float> exp_v <- dist_to_nonkendall_blockgruops_from_my_home collect exp(each* beta_work_allocation);
////				// debug
//				if this_in_kendall {
//					exp_v <- dist_to_nonkendall_blockgruops_from_my_home collect exp(each* (-beta_work_allocation));
//				}
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
	
	action new_construction {
		
	}
	
	action kendall_statistics {
		kendall_occupancy <- kendall_virtual_block.crt_total_pop / sum(landuse collect each.resident_capacity);
		kendall_diversity <- kendall_virtual_block.crt_diversity;
		kendall_low_inc_ratio <- kendall_virtual_block.crt_low_inc_pop / kendall_virtual_block.crt_total_pop;
		kendall_building_energy <- 0.0;
		kendall_density <- kendall_virtual_block.crt_pop_density;
		kendall_finance <- 0.0;
		
		// commute calculation
		mean_commute_distance['all'] <- sum(people collect (each.represented_nb_people * each.commute_distance)) / sum(people collect each.represented_nb_people);
		list<people> people_work_in_kendall <- people where each.work_in_kendall;
		list<people> people_live_in_kendall <- people where each.live_in_kendall;
		list<people> people_work_or_live_in_kendall <- people where (each.work_in_kendall or each.live_in_kendall);
		mean_commute_distance['work_in_kendall'] <- sum(people_work_in_kendall collect (each.represented_nb_people * each.commute_distance)) / 
			sum(people_work_in_kendall collect each.represented_nb_people);
		mean_commute_distance['live_in_kendall'] <- sum(people_live_in_kendall collect (each.represented_nb_people * each.commute_distance)) / 
			sum(people_live_in_kendall collect each.represented_nb_people);
		mean_commute_distance['work_or_live_in_kendall'] <- sum(people_work_or_live_in_kendall collect (each.represented_nb_people * each.commute_distance)) / 
			sum(people_work_or_live_in_kendall collect each.represented_nb_people);
			
		move_max_value_record_list << max(move_to_kendall_pop.values + move_out_of_kendall_pop);  
//		move_max_value_record_list <- [1,2,3,4,5,6,7,8] collect move_max_value_record_list[max(length(move_max_value_record_list)-each, 0)];  //keep only the last 3 elements
//		crt_move_y_axis_upper_range <- int(quantile(move_max_value_record_list, 0.999999999999999))+10;
		if cycle=4 {remove max(move_max_value_record_list) from:move_max_value_record_list;}
		crt_move_y_axis_upper_range <- max(move_max_value_record_list) + 15;
//		write string(crt_move_y_axis_upper_range) + ', ' +  string(int(quantile(move_max_value_record_list, 1))+10) + ', ' + string(max(move_max_value_record_list));
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
	float finance;
	float revene;
	float expenditure;
	init {
	}
}


species blockgroup {
	string geoid;
	float myarea;
	int init_low_inc_pop;
	int init_high_inc_pop;
	float small_size_apt_ratio;
	float rent;
	int crt_nb_available_bedrooms;
	int crt_low_inc_pop;
	int crt_high_inc_pop;
	int crt_total_pop;
	float crt_mean_income;
	float crt_pop_density;
	float crt_diversity;
	list<people> attached_residents <- [];
	map<string,int> nb_workers_in_kendall <- ['low'::0, 'high'::0];
	float pop_density;
	bool is_abstract_whole_kendall <- false;
	
	
	action update_current_population {
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
		if self = kendall_virtual_block { // update kendall virtual block
			crt_nb_available_bedrooms <- sum(landuse collect each.crt_nb_available_bedrooms);
			if crt_nb_available_bedrooms > 0 {
				rent <- sum(landuse collect (each.rent * each.crt_nb_available_bedrooms)) / crt_nb_available_bedrooms;
			} else {
				rent <- 0.0;
			}
			if init_low_inc_pop + init_high_inc_pop > 0 {
				small_size_apt_ratio <- init_low_inc_pop / (init_low_inc_pop + init_high_inc_pop);
			} else {
				small_size_apt_ratio <- 0.0;
			}
			if crt_total_pop > 0{
				crt_diversity <- -(crt_low_inc_pop/crt_total_pop) * log(crt_low_inc_pop/crt_total_pop) - (crt_high_inc_pop/crt_total_pop) * log(crt_high_inc_pop/crt_total_pop);
			} else {
				crt_diversity <- 0.0;
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
						write "Init people warning: a person in kendall block (geoid=" + myself.geoid + ") can not find a grid inside his block with vacant residence, use close available grid instead";
					}
					if home_grid != nil {
						location <- any_location_in(home_grid);
						home_grid.attached_residents << self;
						home_grid.crt_total_pop <- home_grid.crt_total_pop + this_represented_nb_people;
						if this_income = 0 {
							home_grid.crt_low_inc_pop <- home_grid.crt_low_inc_pop + this_represented_nb_people;
						} else {
							home_grid.crt_high_inc_pop <- home_grid.crt_high_inc_pop + this_represented_nb_people;
						}
						home_grid.crt_nb_available_bedrooms <- home_grid.crt_nb_available_bedrooms - this_represented_nb_people;
					} else {
						location <- any_location_in(myself);
						write "Init people warning: all landuse grids are fully occupied, the person has to find a location in nonKendall blocks";
					}
						
				} else {
					location <- any_location_in(myself);
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
	
	action select_residence {
		float t1 <- machine_time;
		list<string> alternative_geoid_list <- [home_block.geoid]; //current blockgroup will always be in the choice set
		bool kendall_in_choice_set <- false;
		if kendall_virtual_block.crt_nb_available_bedrooms > 0 {  // kendall will be in the choice set as long as there are vacant houses 
			if (distance_to_blockgroups_from_my_work_map['250173523001'] < 3000) or (
				distance_to_blockgroups_from_my_work_map['250173523001'] >= 3000 and 
				flip(nb_nonkendall_blocks_in_hom_loc_choice/length(nonkendall_geoid_list+1))
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
		list<float> pop_density_in_these_blocks <- alternative_geoid_list collect (blockgroup_lookup_map[each].crt_pop_density);
		list<float> inc_disparity_to_these_blocks <- alternative_geoid_list collect (abs(blockgroup_lookup_map[each].crt_mean_income - income));
		matrix<float> tmp_x <- matrix([move_or_stay_to_these_blocks, commute_dist_to_these_blocks, size_in_these_blocks,
								rent_in_these_blocks, pop_density_in_these_blocks, inc_disparity_to_these_blocks
		]);
		list<float> v <- list(tmp_x.b_mat);
		list<float> p;
		if represented_nb_people > lower_nb_represented_persons {// if represented_nb_people is big enough, use the real aggreated prob 
			 p <- v_to_p(v);
		} else { //if too few people make choice, aggreated prob may cause round error, like round([0.3,0.3,0.4] * 1 person) = [0,0,0]
//			int chosen_idx <- to_choose(v); // to avoid this round error, use simulated prob by Monte Carlo method
			int chosen_idx <- (v index_of max(v)); //to avoid this round error, always choose the alt with max_v
			p <- v collect 0.0;
			p[chosen_idx] <- 1.0;
		}
		
//////		// debug
//		if work_in_kendall and represented_nb_people > 0 and (live_in_kendall=false) {
//			write '\n------------------';
//			write "Kendall vacant number of residences: " + kendall_virtual_block.crt_nb_available_bedrooms ;
////			write "I am working in Kendall: " + work_in_kendall;
//			write 'Kendall in choice set: ' + kendall_in_choice_set;
//			write "dist: "  + commute_dist_to_these_blocks;
//			write "rent: " + rent_in_these_blocks;
// 			if commute_dist_to_these_blocks[1] != min(commute_dist_to_these_blocks) {
//				write "work grid: " + (work_grid=nil? work_block.geoid : work_grid.id) + ", crt live block: " + home_block.geoid;
//			}
//			write 'My name: ' + self + ', represent number of people: ' + represented_nb_people;
//			if kendall_in_choice_set {write "Kendall prob: " + p[1]; write "Expected move in people: " + string(represented_nb_people*p[1]) + ', round to: ' + round(represented_nb_people*p[1]);}
//		}
//		
		loop idx from:0 to:length(p)-1 {
			int nb_people_to_make_this_choice <- round(represented_nb_people * p[idx]);
			if idx = 0 {
				// choose to stay, pass
			} else if kendall_in_choice_set and idx = 1 and nb_people_to_make_this_choice > 0 { 
				// choose to move to kendall, do a secenod choice on sepecific landuse grid
				list<landuse> alternative_grid_list <- landuse where (each.crt_nb_available_bedrooms > 0);
				if length(alternative_grid_list) = 0 {break;}
				list<float> move_or_stay_to_these_grids <- alternative_grid_list collect 1.0;
				list<float> commute_dist_to_these_grids <- alternative_grid_list collect (distance_to_grids_from_my_work_map[each.id] / 1000);
				list<float> size_in_these_grids <- alternative_grid_list collect (float(each.is_affordable));
				list<float> rent_in_these_grids <- alternative_grid_list collect (each.rent);
				list<float> pop_density_in_these_grids <- alternative_grid_list collect (each.crt_pop_density);
				list<float> inc_disparity_to_these_grids <- alternative_grid_list collect (abs(each.crt_mean_income - income));
				matrix<float> tmp_x_grids <- matrix([move_or_stay_to_these_grids, commute_dist_to_these_grids, size_in_these_grids,
										rent_in_these_grids, pop_density_in_these_grids, inc_disparity_to_these_grids
				]);
				list<float> v_grids <- list(tmp_x_grids.b_mat);
				list<float> p_grids;
				if nb_people_to_make_this_choice > lower_nb_represented_persons {// if nb_people_to_make_this_choice is big enough, use the real aggreated prob 
					 p_grids <- v_to_p(v_grids);
				} else { //if too few people make choice, aggreated prob may cause round error, like round([0.3,0.3,0.4] * 1 person) = [0,0,0]
//					int chosen_idx <- to_choose(v_grids);   // to avoid this round error, use simulated prob by Monte Carlo method
					int chosen_idx <- v_grids index_of max(v_grids);  // to avoid this roud error, always choose the alt with max v
					p_grids <- v_grids collect 0.0;
					p_grids[chosen_idx] <- 1.0;
				}
				
				loop grid_idx from:0 to:length(alternative_grid_list)-1 {
					int nb_people_to_move_to_this_grid <- round(nb_people_to_make_this_choice * p_grids[grid_idx]);
					if nb_people_to_move_to_this_grid > 0 {
						landuse destination_grid <- alternative_grid_list[grid_idx];
						// in case that the actually number of vacant bedrooms is less than the number of people going to move to this grid
						nb_people_to_move_to_this_grid <- min(nb_people_to_move_to_this_grid, destination_grid.crt_nb_available_bedrooms);
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
							each.income=income and each.work_block=work_block and each.work_grid=work_grid
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
					}
				} 
				
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
	
	people copy_myself (blockgroup new_home_block, landuse new_home_grid, point new_home_loc, int new_represented_nb_people, int new_settled_time) {
		create people number: 1 returns: copied_people {
			represented_nb_people <- (new_represented_nb_people = -1) ? myself.represented_nb_people : new_represented_nb_people;
			show_size <- 5 + (20-5)*(represented_nb_people-1)/(50-1);
			mycolor <- myself.mycolor;
			income <- myself.income;
			live_in_kendall <- (new_home_block in kendall_blockgroups) or (new_home_block = 'kendall');
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
			if work_grid != nil and home_grid != nil {commute_distance <- grid_distance_matrix_map[home_grid.id][work_grid.id];}
			else if work_grid = nil and home_grid = nil {commute_distance <- blockgroup_distance_matrix_map[home_block.geoid][work_block.geoid];}
			else if work_grid != nil {commute_distance <- blockgroup_to_grid_distance_matrix_map[home_block.geoid][work_grid.id];}
			else if home_grid != nil {commute_distance <- blockgroup_to_grid_distance_matrix_map[work_block.geoid][home_grid.id];}
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
	float rent_discount_ratio <- 1.0;
	float rent_shift <- 0.0;
	float rent update:max(rent_base*rent_discount_ratio + rent_shift, 0.0);
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
	int nb_workers <- 0;
	list<people> attached_workers;
	bool is_affordable;
	blockgroup associated_blockgroup;
	list<landuse> my_grid_neighors;
	
	init {
		if usage = 'R' {
			is_affordable <- scale='S';
			resident_capacity <- int(myarea * far / (is_affordable ? 40 : 80));
			mycolor <- is_affordable ? #red : #orange;
		}
		if usage = 'O' {
			worker_capacity <- int(myarea * far / (scale='S' ? 20 : 40));
			mycolor <- #blue;
		}
		if usage = 'vacant' {
			mycolor <- #grey;
		}
		my_grid_neighors <- landuse where (each.id in grid_neighbors_map[id]);
		crt_nb_available_bedrooms <- resident_capacity;
	}
	
	action get_associated_blockgroup {
		associated_blockgroup <- first(blockgroup where (each covers self));
		if associated_blockgroup = nil {
			associated_blockgroup <- first(blockgroup where (each covers centroid(self)));
		}
	}
	
	action update_current_population {
		crt_low_inc_pop <- sum((attached_residents where (each.income=0)) collect each.represented_nb_people);
		crt_high_inc_pop <- sum((attached_residents where (each.income=1)) collect each.represented_nb_people);
		crt_total_pop <- crt_low_inc_pop + crt_high_inc_pop;
		crt_nb_available_bedrooms <- max(resident_capacity - crt_total_pop, 0);
		int total_pop_around <- sum(my_grid_neighors collect each.crt_total_pop);
		int total_high_inc_pop_around <- sum(my_grid_neighors collect each.crt_high_inc_pop);
		crt_pop_density <- total_pop_around / (0.64 * length(my_grid_neighors));
		if total_pop_around > 0 {
			crt_mean_income <- total_high_inc_pop_around / total_pop_around;
		} else {
			crt_mean_income <- 0.0;
		}
	}
	
	action change_land_use{
		
	}
	
	action change_far {
		
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
//	float minimum_cycle_duration <- 1#sec;
	parameter "Show Blockgroup" var:show_blockgroup category:Viz;
	parameter "Focus on Kendall Square" var:focus_on_kendall category:Viz;
	parameter "Show Legend" var:show_legend category:Viz;
	parameter "Show People" var:show_people among: ['Show all people', 'Show Kendall people only', 'Hide all people'] category:Viz;
	parameter "Show Workplace" var:show_workplace category:Viz;
	parameter "Show Commute" var:show_commute category:Viz;
	parameter "Show Move" var:show_move category:Viz;
	
	
	parameter "Monthly Subsidy for Low Income People" var:subsidy_low_inc category:"Incentive Policy";
	parameter "Lump-Sum Subsidy for Commute Deccreasing Move" var:sbusidy_less_commute category:"Incentive Policy";
	
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
			chart "Moving Population" type: series background: #white position: {0,0.0} size: {0.5,0.5} y_range:{0,crt_move_y_axis_upper_range}{
				data "Move to Kendall: Low Income" value: move_to_kendall_pop[0] color: #green;
				data "Move to Kendall: High Income" value: move_to_kendall_pop[1] color: #blue;
				data "Move out of Kendall: Low Income" value:move_out_of_kendall_pop[0] color: #red;
				data "Move out of Kendall: High Income" value:move_out_of_kendall_pop[1] color: #orange;
			}
			chart "Occupancy" type:series background: #white position:{0,0.5} size:{0.5,0.5}{
				data "Occupancy" value: kendall_occupancy color: #red;
			}
			chart "Commute Distance" type:series background: #white position:{0.5,0} size:{0.5,0.5} {
				data "All" value: mean_commute_distance['all'] color:#red;
				data "Work in Kendall" value: mean_commute_distance['work_in_kendall'] color:#green;
				data "Live in Kendall" value: mean_commute_distance['live_in_kendall'] color:#blue;
				data "Work or Live in Kendall" value: mean_commute_distance['work_or_live_in_kendall'] color:#yellow;
			}
			chart "Income Structure" type:series background: #white position:{0.5,0.5} size:{0.5,0.5}{
				data "Diversity" value: kendall_diversity color: #red;
				data "Low Income Proportion" value: kendall_low_inc_ratio color: #blue;
			}
		}
//		monitor "Total population" value: sum(people collect (each.represented_nb_people));
//		monitor "Total population moved" value: sum((people where (each.settled_time<1)) collect (each.represented_nb_people));
		inspect name:'Kendall' value:kendall_virtual_block type:agent;
	}
	
}