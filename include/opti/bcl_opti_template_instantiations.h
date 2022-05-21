// (c) Copyright BCL @ Vanderbilt University 2014
// (c) BCL Homepage: http://www.meilerlab.org/bclcommons
// (c) BCL Code Repository: https://github.com/BCLCommons/bcl
// (c)
// (c) The BioChemical Library (BCL) was originally developed by contributing members of the Meiler Lab @ Vanderbilt University.
// (c)
// (c) The BCL is now made available as an open-source software package distributed under the permissive MIT license,
// (c) developed and maintained by the Meiler Lab at Vanderbilt University and contributing members of the BCL Commons.
// (c)
// (c) External code contributions to the BCL are welcome. Please visit the BCL Commons GitHub page for information on how you can contribute.
// (c)
// (c) This file is part of the BCL software suite and is made available under the MIT license.
// (c)

#ifndef BCL_OPTI_TEMPLATE_INSTANTIATIONS_H_
#define BCL_OPTI_TEMPLATE_INSTANTIATIONS_H_

// include the namespace header
#include "bcl_opti.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_opti_criterion_all.h"
#include "bcl_opti_criterion_combine.h"
#include "bcl_opti_criterion_number_iterations.h"
#include "bcl_opti_criterion_unimproved.h"
#include "bcl_opti_ensemble_node.h"
#include "bcl_opti_optimization_identity.h"
#include "bcl_opti_pipeline.h"
#include "assemble/bcl_assemble_protein_model.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace opti
  {

    // Note: If desired (e.g. for purposes of dynamic selection on the command line), a list of all
    // Criterions can be created with this command:
    // grep 'i::Criterion[A-Za-z \n]*<[ \n]*[a-z][^_]' ./source/*/* ./include/*/* ./apps/* ./apps/*/* ./apps/*/*/*
    //      ./apps/*/*/*/* ./example/* ./example/*/* | sed 's/^.*\.[chp]*:[ \t]*//' | sed 's/\([^-]\)>[^>]*$/\1>/' |
    //      sed 's/new //' | sed 's/.*opti::/opti::/' | sed 's/> >.*/>/' | sort | uniq | sed 's/opti:://' |
    //      awk '{print "     BCL_EXPIMP_TEMPLATE template class BCL_API "$0";"}'

    BCL_EXPIMP_TEMPLATE template class BCL_API CriterionAll< assemble::ProteinModel, double>;

    BCL_EXPIMP_TEMPLATE template class BCL_API CriterionCombine< assemble::ProteinModel, double>;

    BCL_EXPIMP_TEMPLATE template class BCL_API CriterionNumberIterations< assemble::ProteinModel, double>;

    BCL_EXPIMP_TEMPLATE template class BCL_API CriterionUnimproved< assemble::ProteinModel, double>;

    BCL_EXPIMP_TEMPLATE template class BCL_API EnsembleNode< assemble::ProteinModel>;

    BCL_EXPIMP_TEMPLATE template class BCL_API OptimizationIdentity< assemble::ProteinModel>;

    BCL_EXPIMP_TEMPLATE template class BCL_API OptimizationIdentity< assemble::Ensemble< assemble::ProteinModel> >;

    BCL_EXPIMP_TEMPLATE template class BCL_API Pipeline< assemble::ProteinModel>;

  } // namespace opti
} // namespace bcl

#endif // BCL_OPTI_TEMPLATE_INSTANTIATIONS_H_
