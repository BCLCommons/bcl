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

// include example header
#include "example.h"
// include the header of the class which this example is for
#include "chemistry/bcl_chemistry_bond_lengths.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_fragment_complete.h"
#include "io/bcl_io_file.h"
#include "linal/bcl_linal_vector_3d_operations.h"
#include "sdf/bcl_sdf_fragment_factory.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_chemistry_bond_lengths.cpp
  //! @details Tests ChemistryBondLengths class which contains small molecule configuration data
  //!
  //! @author mendenjl
  //! @date July 9, 2012
  //! @remarks status complete
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleChemistryBondLengths :
    public ExampleInterface
  {
  public:

    ExampleChemistryBondLengths *Clone() const
    {
      return new ExampleChemistryBondLengths( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    int Run() const
    {

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // construct from conformation interface
      io::IFStream input_sdf;
      const std::string taxol_filename( AddExampleInputPathToFilename( e_Chemistry, "taxol.sdf"));
      BCL_ExampleMustOpenInputFile( input_sdf, taxol_filename);

      // load information into taxol
      chemistry::FragmentComplete taxol;
      taxol = sdf::FragmentFactory::MakeFragment( input_sdf);
      // close the input stream
      io::File::CloseClearFStream( input_sdf);

    /////////////////
    // data access //
    /////////////////

      // Compute the max deviation between the actual bonds lengths and the average bond length known by this
      // class
      double max_deviation( 0.0);
      double max_deviation_outside_ring_oh( 0.0);
      for
      (
        iterate::Generic< const chemistry::AtomConformationalInterface> itr( taxol.GetAtomsIterator());
        itr.NotAtEnd();
        ++itr
      )
      {
        const linal::Vector3D &position( itr->GetPosition());
        for
        (
          storage::Vector< chemistry::BondConformational>::const_iterator
            itr_bond( itr->GetBonds().Begin()), itr_bond_end( itr->GetBonds().End());
          itr_bond != itr_bond_end;
          ++itr_bond
        )
        {
          const linal::Vector3D &position_b( itr_bond->GetTargetAtom().GetPosition());
          const double actual_distance( linal::Distance( position, position_b));
          const double expected_distance
          (
            chemistry::BondLengths().GetBondLength
            (
              itr->GetAtomType(),
              itr_bond->GetBondType()->GetBondData( chemistry::ConfigurationalBondTypeData::e_BondOrderOrAromatic),
              itr_bond->GetTargetAtom().GetAtomType()
            )
          );
          const double deviation( math::Absolute( actual_distance - expected_distance));
          BCL_MessageVrb
          (
            "Deviation between " + itr->GetAtomType().GetName()
            + " " + itr_bond->GetBondType().GetName()
            + " " + itr_bond->GetTargetAtom().GetAtomType().GetName()
            + " = " + util::Format()( deviation)
          );
          max_deviation = std::max( max_deviation, deviation);

          // compute max deviation for bonds outside rings; also skip oxygens and hydrogens, which both have relatively
          // large deviations in bond lengths
          if( !itr_bond->GetBondType()->IsBondInRing())
          {
            if
            (
              itr->GetElementType()->GetAtomicNumber() != size_t( 1)
              && itr->GetElementType()->GetAtomicNumber() != size_t( 8)
              && itr_bond->GetTargetAtom().GetElementType()->GetAtomicNumber() != size_t( 1)
              && itr_bond->GetTargetAtom().GetElementType()->GetAtomicNumber() != size_t( 8)
            )
            {
              max_deviation_outside_ring_oh = std::max( max_deviation_outside_ring_oh, deviation);
            }
          }
        }
      }

      BCL_ExampleIndirectCheck
      (
        max_deviation < 0.25,
        true,
        "max deviation of actual bond lengths with expected bonds lengths in taxol was "
        + util::Format()( max_deviation)
      );

      BCL_ExampleIndirectCheck
      (
        max_deviation_outside_ring_oh < 0.07,
        true,
        "max deviation of actual bond lengths with expected bonds lengths, outside rings and for non O/H bonds in taxol was "
        + util::Format()( max_deviation_outside_ring_oh)
      );

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleSmallMolecule

  const ExampleClass::EnumType ExampleChemistryBondLengths::s_Instance
  (
    GetExamples().AddEnum( ExampleChemistryBondLengths())
  );

} // namespace bcl
