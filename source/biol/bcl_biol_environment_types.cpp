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
#include "biol/bcl_biol_environment_types.h"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_enumerate.hpp"

// external includes - sorted alphabetically

namespace bcl
{
  namespace biol
  {

    //! @brief construct all EnvironmentTypes
    EnvironmentTypes::EnvironmentTypes() :
      util::Enumerate< EnvironmentTypeData, EnvironmentTypes>( false),
      e_MembraneCore(          AddEnum( "MEMBRANE_CORE",                    EnvironmentTypeData( "MC", "MC", 0, false, double(     10.0)))),
      e_GapCoreTransition(     AddEnum( "MEMBRANE_GAP_CORE_TRANSITION",     EnvironmentTypeData( "G1", "TR", 1, true,  double(      2.5)))),
      e_Transition(            AddEnum( "MEMBRANE_TRANSITION",              EnvironmentTypeData( "TR", "TR", 1, false, double(     10.0)))),
      e_GapTransitionSolution( AddEnum( "MEMBRANE_GAP_TRANSITION_SOLUTION", EnvironmentTypeData( "G2", "SO", 2, true,  double(      2.5)))),
      e_Solution(              AddEnum( "SOLUTION",                         EnvironmentTypeData( "SO", "SO", 2, false, double( 100000.0)))),
      e_MembraneInside(        AddEnum( "MEMBRANE_INSIDE",                  EnvironmentTypeData( "MI", "MC", 0, false, double(     10.0)))),
      e_MembraneOutside(       AddEnum( "MEMBRANE_OUTSIDE",                 EnvironmentTypeData( "ME", "MC", 0, false, double(     10.0)))),
      e_SolutionInside(        AddEnum( "SOLUTION_INSIDE",                  EnvironmentTypeData( "SI", "SO", 2, false, double( 100000.0)))),
      e_SolutionOutside(       AddEnum( "SOLUTION_OUTSIDE",                 EnvironmentTypeData( "SE", "SO", 2, false, double( 100000.0))))
    {
    }

    //! @brief function to deduce EnvironmentType from two letter code
    //! @param TWO_LETTER_CODE two letter code descriptor for environment type of interest
    //! @return EnvironmentType specified by the given TWO_LETTER_CODE
    const EnvironmentType &EnvironmentTypes::EnvironmentTypeFromTwoLetterCode( const std::string &TWO_LETTER_CODE) const
    {
      // iterate over all EnvironmentTypes
      for
      (
        const_iterator type_itr( Begin()), type_itr_end( End());
        type_itr != type_itr_end;
        ++type_itr
      )
      {
        // if ONE_LETTER_CODE matches this EnvironmentTypes's two letter code
        if( ( *type_itr)->GetTwoLetterCode() == TWO_LETTER_CODE)
        {
          // return this EnvironmentType
          return *type_itr;
        }
      }

      // if no match was found return undefined EnvironmentType
      return GetEnvironmentTypes().e_Undefined;
    }

    //! @brief returns number of reduced types
    //! @return number of reduced types
    size_t EnvironmentTypes::GetNumberReducedTypes() const
    {
      // end
      return GetReducedTypes().GetSize();
    }

    //! @brief returns a vector that contains the 3 reduces types, e_MembraneCore, e_Transition and e_Solution
    //! @return a vector that contains the 3 reduces types, e_MembraneCore, e_Transition and e_Solution
    const storage::Vector< EnvironmentType> &EnvironmentTypes::GetReducedTypes() const
    {
      // initialize a static vector to store reduced types
      static const storage::Vector< EnvironmentType> s_reduced_types_vector
      (
        storage::Vector< EnvironmentType>::Create( e_MembraneCore, e_Transition, e_Solution)
      );

      // end
      return s_reduced_types_vector;
    }

    //! @brief construct on access function for all EnvironmentTypes
    //! @return reference to only instance of EnvironmentTypes enum
    const EnvironmentTypes &GetEnvironmentTypes()
    {
      return EnvironmentTypes::GetEnums();
    }

  } // namespace biol

  namespace util
  {

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    template class BCL_API Enumerate< biol::EnvironmentTypeData, biol::EnvironmentTypes>;

  } // namespace util
} // namespace bcl
