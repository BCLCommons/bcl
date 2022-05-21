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
#include "chemistry/bcl_chemistry_conformation_comparison_by_real_space_rmsd.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_conformation_interface.h"
#include "quality/bcl_quality_rmsd.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

    // add the interface to the set of known implementations
    const util::SiPtr< const util::ObjectInterface> ConformationComparisonByRealSpaceRmsd::s_Instance
    (
      util::Enumerated< ConformationComparisonInterface>::AddInstance( new ConformationComparisonByRealSpaceRmsd())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! virtual copy constructor
    ConformationComparisonByRealSpaceRmsd *ConformationComparisonByRealSpaceRmsd::Clone() const
    {
      return new ConformationComparisonByRealSpaceRmsd( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! the class name as const ref std::string
    const std::string &ConformationComparisonByRealSpaceRmsd::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &ConformationComparisonByRealSpaceRmsd::GetAlias() const
    {
      static const std::string s_Name( "RealSpaceRMSD");
      return s_Name;
    }

  /////////////////
  //  operations //
  /////////////////

    //! @brief align two small molecule objects and find the RMSD
    //! @param MOLECULE_A - first molecule being aligned
    //! @param MOLECULE_A - second molecule being aligned
    //! @return the RMSD between first and second molecule
    double ConformationComparisonByRealSpaceRmsd::operator()
    (
      const ConformationInterface &MOLECULE_A,
      const ConformationInterface &MOLECULE_B
    ) const
    {
      // Determine if conformations can be compared
      if( !ConformationComparisonInterface::ConformationsAreComparable( MOLECULE_A, MOLECULE_B))
      {
        // nope, so return an undefined
        return util::GetUndefined< double>();
      }

      // find the RMSDs between the coordinate vectors
      double rmsd( quality::RMSD::RealSpaceRMSD( MOLECULE_A.GetAtomCoordinates(), MOLECULE_B.GetAtomCoordinates()));

      // if the rmsd is very close to zero, set it to zero
      // assume that errors of size ~std::numeric_limits< double>::epsilon() are accumulated in the RMSD calculation for every atom
      if( math::EqualWithinAbsoluteTolerance( 0.0, rmsd, MOLECULE_A.GetNumberAtoms() * std::numeric_limits< double>::epsilon()))
      {
        rmsd = 0.0;
      }

      return rmsd;
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer ConformationComparisonByRealSpaceRmsd::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Calculates the root-mean-squared-deviation between two molecules that are identical on"
        "the constitutional level, and whose atoms in the same order.  No superimposition is performed."
      );
      return parameters;
    }

  } // namespace chemistry
} // namespace bcl
