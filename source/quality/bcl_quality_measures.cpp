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
#include "quality/bcl_quality_measures.h"

// includes from bcl - sorted alphabetically
#include "command/bcl_command_flag_dynamic.h"
#include "command/bcl_command_parameter.h"
#include "command/bcl_command_parameter_check_enumerate.h"
#include "quality/bcl_quality_average.h"
#include "quality/bcl_quality_dme.h"
#include "quality/bcl_quality_dmf.h"
#include "quality/bcl_quality_gdt.h"
#include "quality/bcl_quality_rmsd.h"
#include "quality/bcl_quality_superimpose_measures.h"
#include "util/bcl_util_enumerate.hpp"

// external includes - sorted alphabetically

namespace bcl
{
  namespace quality
  {

  //////////
  // data //
  //////////

    //! @brief returns the default distance cutoff set for HA high accuracy
    //! @return the default distance cutoff set HA
    const storage::Set< double> &Measures::GetDistanceCutoffsHA()
    {
      // initialize default static distance cutoff vector
      static const storage::Set< double> s_cutoffs( storage::Set< double>::Create( 0.5, 1.0, 2.0, 4.0));
      return s_cutoffs;
    }

    //! @brief returns the default distance cutoff set for TS (total score)
    //! @return the default distance cutoff set TS
    const storage::Set< double> &Measures::GetDistanceCutoffsTS()
    {
      // initialize default static distance cutoff vector
      static const storage::Set< double> s_cutoffs( storage::Set< double>::Create( 1.0, 2.0, 4.0, 8.0));
      return s_cutoffs;
    }

    //! @brief return command line flag for defining the quality measures to be calculated
    //! @return command line flag for defining the quality measures to be calculated
    util::ShPtr< command::FlagInterface> &Measures::GetFlagQualityMeasures()
    {
      // initialize static flag
      static util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagDynamic
        (
          "quality",
          "\tlist of quality measures to be calculated",
          command::Parameter
          (
            "quality_measure",
            "\tquality measure to be calculated",
            command::ParameterCheckEnumerate< Measures>()
          ),
          0,
          GetMeasures().GetEnumCount()
        )
      );

      // end
      return s_flag;
    }

    //! @brief function to return the list of quality measures in a set defined by the command line flag
    //! @return the list of quality measures in a set defined by the command line flag
    storage::Set< Measure> Measures::GetCommandLineQualityMeasures()
    {
      return GetFlagQualityMeasures()->GetObjectSet< Measure>();
    }

    //! @brief construct all Measures
    Measures::Measures() :
      e_RMSD(                   AddEnum( "RMSD"                  , *GetSuperimposeMeasures().e_RMSD)),
      e_RMSD_NoSuperimposition( AddEnum( "RMSD_NoSuperimposition", util::ShPtr< MeasureInterface>( new RMSD( false)))),
      e_RMSD_XYSuperimposition( AddEnum( "RMSD_XYSuperimposition", *GetSuperimposeMeasures().e_RMSD_XYSuperimposition)),
      e_DME(                    AddEnum( "DME"                   , util::ShPtr< MeasureInterface>( new DME()))),
      e_DMF_HA(                 AddEnum( "DMF_HA"                , util::ShPtr< MeasureInterface>( new DMF( GetDistanceCutoffsHA())))),
      e_DMF_TS(                 AddEnum( "DMF_TS"                , util::ShPtr< MeasureInterface>( new DMF( GetDistanceCutoffsTS())))),
      e_LCS(                    AddEnum( "LCS"                   , *GetSuperimposeMeasures().e_LCS)),
      e_GDT_HA(                 AddEnum( "GDT_HA"                , util::ShPtr< MeasureInterface>( GDT::CreateAverageGDT( GetDistanceCutoffsHA()).Clone()))),
      e_GDT_TS(                 AddEnum( "GDT_TS"                , util::ShPtr< MeasureInterface>( GDT::CreateAverageGDT( GetDistanceCutoffsTS()).Clone()))),
      e_GDT_1A(                 AddEnum( "GDT_1A"                , *GetSuperimposeMeasures().e_GDT_1A)),
      e_GDT_2A(                 AddEnum( "GDT_2A"                , *GetSuperimposeMeasures().e_GDT_2A)),
      e_GDT_4A(                 AddEnum( "GDT_4A"                , *GetSuperimposeMeasures().e_GDT_4A)),
      e_GDT_8A(                 AddEnum( "GDT_8A"                , *GetSuperimposeMeasures().e_GDT_8A)),
      e_MaxSub(                 AddEnum( "MaxSub"                , *GetSuperimposeMeasures().e_MaxSub)),
      e_Zero(                   AddEnum( "Zero"                  , *GetSuperimposeMeasures().e_NoSuperimpose))
    {
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &Measures::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief function that returns the static instance of the Measures class
    //! @return the static instance of the Measures class
    Measures &GetMeasures()
    {
      return Measures::GetEnums();
    }

  } // namespace quality

  namespace util
  {

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    template class BCL_API Enumerate< ShPtr< quality::MeasureInterface>, quality::Measures>;

  } // namespace util
} // namespace bcl
