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

// initialize the static initialization fiasco finder, if macro ENABLE_FIASCO_FINDER is defined
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "quality/bcl_quality_superimpose_measures.h"

// includes from bcl - sorted alphabetically
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter.h"
#include "command/bcl_command_parameter_check_enumerate.h"
#include "quality/bcl_quality_const_measure.h"
#include "quality/bcl_quality_gdt.h"
#include "quality/bcl_quality_lcs.h"
#include "quality/bcl_quality_maxsub.h"
#include "quality/bcl_quality_rmsd.h"
#include "util/bcl_util_enumerate.hpp"

// external includes - sorted alphabetically

namespace bcl
{
  namespace quality
  {
  //////////
  // data //
  //////////

    //! @brief return command line flag for defining the measures to calculate superimposition
    //! @return command line flag for defining the measures to calculate superimposition
    util::ShPtr< command::FlagInterface> &SuperimposeMeasures::GetFlagSuperimposeMeasure()
    {
      // initialize static flag
      static util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagStatic
        (
          "superimpose",
          "\tFlag for defining the quality measure to use for superimposing models onto template/native model",
          command::Parameter
          (
            "superimpose_measure",
            "\tquality measure to use for superimposition",
            command::ParameterCheckEnumerate< SuperimposeMeasures>(),
            GetSuperimposeMeasures().e_NoSuperimpose.GetName()
          )
        )
      );

      // end
      return s_flag;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief construct all Measures
    SuperimposeMeasures::SuperimposeMeasures() :
      e_RMSD(                   AddEnum( "RMSD",                   util::ShPtr< SuperimposeInterface>( new RMSD( true)))),
      e_RMSD_XYSuperimposition( AddEnum( "RMSD_XYSuperimposition", util::ShPtr< SuperimposeInterface>( new RMSD( true, true)))),
      e_GDT_1A(                 AddEnum( "GDT_1A",                 util::ShPtr< SuperimposeInterface>( new GDT( 1.0)))),
      e_GDT_2A(                 AddEnum( "GDT_2A",                 util::ShPtr< SuperimposeInterface>( new GDT( 2.0)))),
      e_GDT_4A(                 AddEnum( "GDT_4A",                 util::ShPtr< SuperimposeInterface>( new GDT( 4.0)))),
      e_GDT_8A(                 AddEnum( "GDT_8A",                 util::ShPtr< SuperimposeInterface>( new GDT( 8.0)))),
      e_LCS(                    AddEnum( "LCS",                    util::ShPtr< SuperimposeInterface>( new LCS()))),
      e_MaxSub(                 AddEnum( "MaxSub",                 util::ShPtr< SuperimposeInterface>( new MaxSub()))),
      e_NoSuperimpose(          AddEnum( "NoSuperimpose",          util::ShPtr< SuperimposeInterface>( new ConstMeasure())))
    {
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &SuperimposeMeasures::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief function that returns static instance of SuperimposeMeasures
    //! @return static instance of SuperimposeMeasures
    SuperimposeMeasures &GetSuperimposeMeasures()
    {
      return SuperimposeMeasures::GetEnums();
    }

  } // namespace quality

  namespace util
  {

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    template class BCL_API Enumerate< ShPtr< quality::SuperimposeInterface>, quality::SuperimposeMeasures>;

  } // namespace util
} // namespace bcl
