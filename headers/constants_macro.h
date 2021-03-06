#pragma once


/////////////////////////////////////
// propgators
/////////////////////////////////////

// 24ID
//


std::string Kaon_four_point_24ID(int traj) {
  // std::string path = "/projects/CSC249ADSE03/yidizhao/KGG_config/24ID/typeI/old/KGG_typeI." + std::to_string(traj);
  // std::string path = "/home/yidizhao/cooley/decays/kaon/KGG_typeI." + std::to_string(traj);
  // std::string path = "/home/yidizhao/cooley/decays/kaon/old_config/KGG_typeI." + std::to_string(traj);
  std::string path = "/hpcgpfs01/work/lqcd/qcdqedta/ydzhao/24ID/KGG/KGG_typeI." + std::to_string(traj);
  assert(dirExists(path));
  return path;
}

std::string point_path_24ID(int traj) {
  std::string path = "/home/ljin/application/Public/Muon-GM2-cc/jobs/24D/discon-1/results/prop-hvp ; results=" + std::to_string(traj) + "/huge-data/prop-point-src";
  assert(dirExists(path));
  return path;
}

std::string point_path_strange_24ID(int traj) {
  std::string path = "/home/ljin/application/Public/Muon-GM2-cc/jobs/24D/discon-1/results/prop-hvp ; results=" + std::to_string(traj) + "/huge-data/prop-point-src";
  assert(dirExists(path));
  return path;
}

std::string wall_path_s_24ID(int traj, int t) {
  std::string path = "/home/ljin/application/Public/Qlat-CPS-cc/jobs/wall-src-strange/results/24D-0.00107/results="+ std::to_string(traj) + "/huge-data/wall_src_propagator/strange ; t=" + std::to_string(t);
  assert(dirExists(path));
  return path;
}

std::string wall_path_l_24ID(int traj, int t) {
  std::string path = "/home/ljin/application/Public/Qlat-CPS-cc/jobs/24D/wall-src/results/results=" + std::to_string(traj) + "/huge-data/wall_src_propagator/t=" + std::to_string(t);
  assert(dirExists(path));
  return path;
}

std::string gauge_transform_path_24D(int traj) {
	return "/home/ljin/application/Public/Qlat-CPS-cc/jobs/24D/wall-src/results/results=" + std::to_string(traj) + "/huge-data/gauge-transform";
}

std::string loop_path_24ID(int traj) {
  std::string path ="/home/ljin/application/Public/Muon-GM2-cc/jobs/convert/data/24D/light-minus-heavy-1/traj=0" + std::to_string(traj); 
  assert(dirExists(path));
  return path;
}

// 32ID

std::string wall_path_strange_32ID(int traj, int t) {
  std::string path = "/home/ljin/application/Public/Qlat-CPS-cc/jobs/wall-src-strange/results/32D-0.00107/results="+ std::to_string(traj) + "/huge-data/wall_src_propagator/strange ; t=" + std::to_string(t);
  assert(dirExists(path));
  return path;
}

std::string wall_path_ud_32ID(int traj, int t) {
  std::string path = "/home/ljin/application/Public/Qlat-CPS-cc/jobs/32D/wall-src/results/32D-0.00107/results=" + std::to_string(traj) + "/huge-data/wall_src_propagator/t=" + std::to_string(t);
  assert(dirExists(path));
  return path;
}


std::string gauge_transform_path_32ID(int traj) {
	return "/home/ljin/application/Public/Qlat-CPS-cc/jobs/32D/wall-src/results/32D-0.00107/results=" + std::to_string(traj) + "/huge-data/gauge-transform";
}


std::string point_path_32ID(int traj) {
	std::string path = "/home/ljin/application/Public/Muon-GM2-cc/jobs/32D/discon-1/results/prop-hvp ; results=" + std::to_string(traj) + "/huge-data/prop-point-src";
  assert(dirExists(path));
  return path;
}

std::string gauge_transform_path(int traj) {
  return "/home/ljin/application/Public/Qlat-CPS-cc/jobs/32D/wall-src/results/32D-0.00107/results=" + std::to_string(traj) + "/huge-data/gauge-transform";
}


// 48I


std::string wall_path_ud_48I(int traj, int t) {
   std::string path = "/home/ljin/application/Public/Qlat-CPS-cc/jobs/48I/wall-src/results/48I-0.00078/results=" + std::to_string(traj) + "/huge-data/wall_src_propagator/t=" + std::to_string(t);
  assert(dirExists(path));
  return path;
}


/////////////////////////////////////
// Pion three point functions
////////////////////////////////////


std::string three_point_path(int traj, const std::string &ensemble, const std::string &target) {
  std::string path;
  if(ensemble=="Pion_24ID") path = "/home/ljin/application/Public/Qlat-CPS-cc/jobs/em-corr/results/24D-0.00107/results=" + std::to_string(traj) + "/contraction-with-point/pion_gg/";    
  else if(ensemble=="Pion_32ID") path = "/home/ljin/application/Public/Qlat-CPS-cc/jobs/em-corr/results/32D-0.00107/results=" + std::to_string(traj) + "/contraction-with-point/pion_gg/";    
  else if(ensemble=="Pion_32IDF") path = "/home/ljin/application/Public/Qlat-CPS-cc/jobs/em-corr/results/32Dfine-0.0001/results=" + std::to_string(traj) + "/contraction-with-point/pion_gg/";    
  else if(ensemble=="Pion_48I") path = "/home/ljin/application/Public/Qlat-CPS-cc/jobs/em-corr/results/48I-0.00078/results=" + std::to_string(traj) + "/contraction-with-point/pion_gg/";   // the old three point functions that Luchang accidentally deleted 
  else if(ensemble=="Pion_48I_pqpm") path = "/home/ljin/application/Public/Muon-GM2-cc/jobs/final-run/48I/run/results/48I-0.00078/results=" + std::to_string(traj) + "/contraction-with-points/pion_gg/";    
  else if(ensemble=="Pion_64I") path = "/home/ljin/application/Public/Muon-GM2-cc/jobs/final-run/64I/run/results/64I-0.000678/results="  + std::to_string(traj) +  "/contraction-with-points/pion_gg/";    
  else throw("ensemble unknown");

  if(target=="decay" || target=="fission" || target=="decay_cheng") path += target;
  else throw("target unknown");

  if(ensemble!="Pion_64I" && ensemble!="Pion_48I_pqpm") assert(dirExists(path)); // for 64I and 48I_pqpm, I have to add "_type_1" and "_type_2" to the file name later.
  return path;
}


std::string three_point_disc_24ID(int traj) {
  std::string path = "/home/yidizhao/cooley/pGG_config/24ID/disc/t-min=20/pGG_disc." + std::to_string(traj);
  // std::string path = "/projects/CSC249ADSE03/yidizhao/pGG_config/24ID/disc/t-min=20/pGG_disc." + std::to_string(traj);
  // std::string path = "/projects/CSC249ADSE03/yidizhao/pGG_config/24ID/disc/t-min=10/pGG_disc." + std::to_string(traj);
  assert(dirExists(path));
  return path;
}

std::string three_point_disc2_32ID(int traj) {
  std::string path = "/gpfs/mira-fs0/projects/CSC249ADSE03/yidizhao/pGG_config/32ID/disc_2/pGG_disc." + std::to_string(traj);
  assert(dirExists(path));
  return path;
}

// std::string three_point_24ID(int traj) {
//   // std::string path = "/gpfs/mira-fs1/projects/HadronicLight_4/ctu/hlbl/hlbl-pion/TwoPointWallCorrField/24D-0.00107/ama/t-min=0010/results=" + std::to_string(traj) + "/avg ; type=0";    
//   std::string path = "/home/ljin/application/Public/Qlat-CPS-cc/jobs/em-corr/results/24D-0.00107/results=" + std::to_string(traj) + "/contraction-with-point/pion_gg/decay_cheng";    
//   assert(dirExists(path));
//   return path;
// }
//
// std::string three_point_32ID(int traj) {
// 	// std::string str_traj = std::to_string(traj);
// 	// if(str_traj.size()==3) str_traj = "0" + str_traj;
//   //
//   // std::string path = "/gpfs/mira-fs1/projects/HadronicLight_4/ctu/hlbl/hlbl-pion/TwoPointWallCorrField/32D-0.00107/ama/t-min=0010/results=" + str_traj + "/avg ; type=0";    
//
//   std::string path = "/home/ljin/application/Public/Qlat-CPS-cc/jobs/em-corr/results/32D-0.00107/results=" + std::to_string(traj) + "/contraction-with-point/pion_gg/decay_cheng";    
//   // std::string path = "/home/ljin/application/Public/Qlat-CPS-cc/jobs/em-corr/results/32D-0.00107/results=" + std::to_string(traj) + "/contraction-with-point/pion_gg/decay";    
//   assert(dirExists(path));
//   return path;
// }
//
//
// std::string three_point_path_32IDF(int traj) {
//
//   std::string path = "/home/ljin/application/Public/Qlat-CPS-cc/jobs/em-corr/results/32Dfine-0.0001/results="+ std::to_string(traj) + "/contraction-with-point/pion_gg/decay_cheng";
//   assert(dirExists(path));
// 	return path;
// }
//
// std::string three_point_48I(int traj) {
//
//   std::string path = "/home/ljin/application/Public/Qlat-CPS-cc/jobs/em-corr/results/48I-0.00078/results=" + std::to_string(traj) + "/contraction-with-point/pion_gg/decay_cheng";
//   assert(dirExists(path));
// 	return path;
// }
//



