


int mpsxx::qcdmrg::Rundmrg (mpsxx::qcdmrg::DMRG_input& input) 
{
  if(!input.restart()) mpsxx::qcdmrg::Build_random_system(input);

  for(auto id = input.schedule().begin(); id != input.schedule().end(); ++id)
  {
    mpsxx::qcdmrg::Sweep_optimization(*id);
  }
}
