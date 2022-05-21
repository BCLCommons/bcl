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
#include "chemistry/bcl_chemistry_conformation_comparison_by_dihedrals.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_priority_dihedral_angles.h"
#include "io/bcl_io_serialization.h"
#include "math/bcl_math_statistics.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

    // add the interface to the set of known implementations
    const util::SiPtr< const util::ObjectInterface> ConformationComparisonByDihedrals::s_Instance
    (
      util::Enumerated< ConformationComparisonInterface>::AddInstance( new ConformationComparisonByDihedrals())
    );

  ///////////
  // enums //
  ///////////

    //! @brief comparison method as string
    //! @param METHOD the method desired
    //! @return the method as string
    const std::string &ConformationComparisonByDihedrals::GetMethodName( const Method &METHOD)
    {
      static const std::string s_Methods[] =
      {
        "Max",
        "Sum",
        "Mean",
        "Norm",
        GetStaticClassName< Method>()
      };
      return s_Methods[ METHOD];
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    ConformationComparisonByDihedrals::ConformationComparisonByDihedrals( const Method &METHOD) :
      m_Method( METHOD)
    {
    }

    //! virtual copy constructor
    ConformationComparisonByDihedrals *ConformationComparisonByDihedrals::Clone() const
    {
      return new ConformationComparisonByDihedrals( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! the class name as const ref std::string
    const std::string &ConformationComparisonByDihedrals::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &ConformationComparisonByDihedrals::GetAlias() const
    {
      static const std::string s_Name( "Dihedral");
      return s_Name;
    }

  /////////////////
  //  operations //
  /////////////////

    //! @brief align two small molecule objects and find the RMSD
    //! @param MOLECULE_A - first molecule being aligned
    //! @param MOLECULE_A - second molecule being aligned
    //! @return the RMSD between first and second molecule
    double ConformationComparisonByDihedrals::operator()
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

      // if there were 0-2 atoms in the molecules, then there are no dihedral angles, so the difference is always 0
      if( MOLECULE_A.GetNumberAtoms() <= size_t( 3))
      {
        return 0.0;
      }
      // calculate the absolute difference in the dihedral angles from the conformations of each molecule
      linal::Vector< double> abs_diff_dihedrals( PriorityDihedralAngles()( MOLECULE_A, MOLECULE_B));

      // if there are no dihedral angles, return 0
      if( abs_diff_dihedrals.IsEmpty())
      {
        return 0.0;
      }

      double result( 0.0);

      // now compute the results of the operation chosen by the user
      switch( m_Method)
      {
        case e_Max: //! Choose the maximum (absolute) difference in dihedral angles
          result = math::Statistics::MaximumValue( abs_diff_dihedrals.Begin(), abs_diff_dihedrals.End());
          break;
        case e_Sum: //! Add the absolute differences in dihedral angles
          result = math::Statistics::Sum( abs_diff_dihedrals.Begin(), abs_diff_dihedrals.End());
          break;
        case e_Mean: //! Get the average absolute difference in dihedral angles
          result = math::Statistics::Mean( abs_diff_dihedrals.Begin(), abs_diff_dihedrals.End());
          break;
        case e_Norm: //! Get the euclidean norm difference in dihedral angles,
          result = math::Statistics::Norm( abs_diff_dihedrals.Begin(), abs_diff_dihedrals.End());
          break;
        default:
          BCL_MessageCrt( " No valid method given for comparison of dihedral angles");
      };

      // if the result is very close to zero, set it to zero
      // assume that errors of size ~std::numeric_limits< double>::epsilon() are accumulated in the calculation for every atom
      if( math::EqualWithinAbsoluteTolerance( 0.0, result, MOLECULE_A.GetNumberAtoms() * std::numeric_limits< double>::epsilon()))
      {
        result = 0.0;
      }

      return result;
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer ConformationComparisonByDihedrals::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Calculates a function of the differences between all dihedral angles in the molecule"
      );
      parameters.AddInitializer
      (
        "method",
        "statistic to calculate on the absolute differences in dihedral angles",
        io::Serialization::GetAgent( &m_Method),
        GetMethodName( e_Max)
      );
      return parameters;
    }

  } // namespace chemistry
} // namespace bcl
