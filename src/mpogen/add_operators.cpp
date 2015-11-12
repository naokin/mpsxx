#include "add_operators.h"

const double mpsxx::OpComponent<mpsxx::op::QC::I>::value[4] = { 1.0, 1.0, 1.0, 1.0 };

const double mpsxx::OpComponent<mpsxx::op::QC::CREA>::value[2] = { 1.0,-1.0 };
const double mpsxx::OpComponent<mpsxx::op::QC::CREB>::value[2] = { 1.0, 1.0 };
const double mpsxx::OpComponent<mpsxx::op::QC::DESA>::value[2] = { 1.0,-1.0 };
const double mpsxx::OpComponent<mpsxx::op::QC::DESB>::value[2] = { 1.0, 1.0 };

const double mpsxx::OpComponent<mpsxx::op::QC::CREA_CREB>::value[1] = {-1.0 };
const double mpsxx::OpComponent<mpsxx::op::QC::CREB_CREA>::value[1] = { 1.0 };
const double mpsxx::OpComponent<mpsxx::op::QC::DESA_DESB>::value[1] = { 1.0 };
const double mpsxx::OpComponent<mpsxx::op::QC::DESB_DESA>::value[1] = {-1.0 };

const double mpsxx::OpComponent<mpsxx::op::QC::CREA_DESA>::value[2] = { 1.0, 1.0 };
const double mpsxx::OpComponent<mpsxx::op::QC::CREB_DESB>::value[2] = { 1.0, 1.0 };

const double mpsxx::OpComponent<mpsxx::op::QC::CREA_DESB>::value[1] = { 1.0 };
const double mpsxx::OpComponent<mpsxx::op::QC::CREB_DESA>::value[1] = { 1.0 };

const double mpsxx::OpComponent<mpsxx::op::QC::CREA_CREB_DESB>::value[1] = {-1.0 };
const double mpsxx::OpComponent<mpsxx::op::QC::CREB_CREA_DESA>::value[1] = { 1.0 };
const double mpsxx::OpComponent<mpsxx::op::QC::CREB_DESB_DESA>::value[1] = {-1.0 };
const double mpsxx::OpComponent<mpsxx::op::QC::CREA_DESA_DESB>::value[1] = { 1.0 };

const double mpsxx::OpComponent<mpsxx::op::QC::CREA_DESA_CREB_DESB>::value[1] = { 1.0 };

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
