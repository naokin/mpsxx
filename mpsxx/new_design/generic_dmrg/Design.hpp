
//
//  What's necessary?
//
//  + Reference configuration space: single det for each q#
//
//  + Configuration space to be projected: small # of dets -> defines "d", i.e. physical index
//


struct Slater {

   /// occupied spin orbitals (ref.)
   std::vector<int> occ_;

   /// virtual spin orbitals (ref.)
   std::vector<int> vir_;

};

std::vector<Slater> make_dets(const FermiQ& _q, const Slater& _Ref)
{
   FermiQ _qRef(_Ref);
   FermiQ _qDif = _q - _qRef;

   long _pDif = _qDif.p();

   if(_pDif == 0) {
   }
   else {
      // Creates IP dets
      // Creates EA dets
   }
}
