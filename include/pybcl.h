// (c) Copyright BCL @ Vanderbilt University 2014
// (c) BCL Homepage: http://www.meilerlab.org/bclcommons
// (c) The BCL software is developed by the contributing members of the BCL @ Vanderbilt University
// (c) This file is part of the BCL software suite and is made available under license.
// (c) To view or modify this file, you must enter into one of the following agreements if you have not done so already:
// (c) For academic and non-profit users:
// (c)   the BCL Academic Single-User License, available at http://www.meilerlab.org/bclcommons/license
// (c) For commercial users:
// (c)   The BCL Commercial Site License, available upon request from bcl-support-commercial@meilerlab.org
// (c) For BCL developers at Vanderbilt University:
// (c)   The BCL Developer Agreement, available at http://www.meilerlab.org/bclcommons/developer_agreement
// (c)
// (c)   As part of all such agreements, this copyright notice must appear, verbatim and without addition, at the
// (c) top of all source files of the BCL project and may not be modified by any party except the BCL developers at
// (c) Vanderbilt University.
// (c)   The BCL copyright and license yields to non-BCL copyrights and licenses where indicated by code comments.
// (c)   Questions about this copyright notice or license agreement may be emailed to bcl-support-academic@meilerlab.org
// (c) (for academic users) or bcl-support-commercial@meilerlab.org (for commercial users)

#ifndef PYBCL_H_
#define PYBCL_H_

// various axulity includes for PyBCL
#include <bcl/include/descriptor/bcl_descriptor_molecule_2da_smooth_sign_code.h>
#include <bcl/include/linal/bcl_linal.h>
#include <bcl/include/util/bcl_util_sh_ptr.h>
#include <bcl/include/util/bcl_util_si_ptr.h>
#include <bcl/include/util/bcl_util_object_interface.h>

#include <bcl/include/cluster/bcl_cluster_output_interface.h>
#include <bcl/include/cluster/bcl_cluster_output_classes.h>
#include <bcl/include/util/bcl_util_enum.h>
#include <bcl/include/util/bcl_util_enumerate.h>
#include <bcl/include/util/bcl_util_enumerate.hpp>
#include <bcl/include/linal/bcl_linal_vector.h>

template class bcl::util::Enumerate<bcl::util::ShPtr<bcl::cluster::OutputInterface<bcl::linal::Vector<float>, float> >, bcl::cluster::OutputClasses<bcl::linal::Vector<float>, float> >;

namespace bcl
{
} // namespace bcl

#endif // PYBCL_H_
