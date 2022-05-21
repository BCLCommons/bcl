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
#include "score/bcl_score_restraint_noe_knowledge_based.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "io/bcl_io_file.h"
#include "math/bcl_math_gnuplot_heatmap.h"
#include "restraint/bcl_restraint_atom_distance.h"
#include "restraint/bcl_restraint_locator_coordinates_hydrogen.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_score_restraint_noe_knowledge_based.cpp
  //! @details scores noe restraints based on a knowledge based potential comparing results for different bond
  //!           numbers.
  //! @author akinlr, weinerbe
  //! @date 01/05/2011
  //! @remarks status complete
  //! @remarks reviewed by weinerbe on 01/29/2011
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleScoreRestraintNoeKnowledgeBased :
    public ExampleInterface
  {
  public:

    ExampleScoreRestraintNoeKnowledgeBased *Clone() const
    { return new ExampleScoreRestraintNoeKnowledgeBased( *this);}

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    virtual const std::string &GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    int Run() const
    {
      // get a protein model
      assemble::ProteinModel protein_model( Proteins::GetModel( AddExampleInputPathToFilename( e_Biology, "1ubi.pdb")));

      // create locators to find atoms to be used in restraints
      util::ShPtr< assemble::LocatorAtomCoordinatesInterface> locator_43_hd
      (
        new restraint::LocatorCoordinatesHydrogen( 'A', 43, biol::GetAtomTypes().HD21)
      );
      util::ShPtr< assemble::LocatorAtomCoordinatesInterface> locator_30_hd
      (
        new restraint::LocatorCoordinatesHydrogen( 'A', 30, biol::GetAtomTypes().HD13)
      );
      util::ShPtr< assemble::LocatorAtomCoordinatesInterface> locator_17_hg
      (
        new restraint::LocatorCoordinatesHydrogen( 'A', 17, biol::GetAtomTypes().HG11)
      );
      util::ShPtr< assemble::LocatorAtomCoordinatesInterface> locator_13_hg
      (
        new restraint::LocatorCoordinatesHydrogen( 'A', 13, biol::GetAtomTypes().HG13)
      );
      util::ShPtr< assemble::LocatorAtomCoordinatesInterface> locator_1_hb
      (
        new restraint::LocatorCoordinatesHydrogen( 'A', 1, biol::GetAtomTypes().HB3)
      );
      util::ShPtr< assemble::LocatorAtomCoordinatesInterface> locator_71_hb
      (
        new restraint::LocatorCoordinatesHydrogen( 'A', 71, biol::GetAtomTypes().HB3)
      );

      // create restraints
      restraint::AtomDistance restraint_a
      (
        locator_43_hd,
        locator_30_hd,
        util::ShPtr< restraint::Distance>( new restraint::Distance( 3.0, 3.75, 1.8))
      );
      restraint::AtomDistance restraint_b
      (
        locator_17_hg,
        locator_13_hg,
        util::ShPtr< restraint::Distance>( new restraint::Distance( 6.0, 7.0, 1.8))
      );
      restraint::AtomDistance restraint_c
      (
        locator_1_hb,
        locator_71_hb,
        util::ShPtr< restraint::Distance>( new restraint::Distance( 5.0, 6.5, 1.8))
      );

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // test Default constructor
      BCL_MessageStd( "test Default def_constructor");
      score::RestraintNoeKnowledgeBased def_constr;
      BCL_Example_Check
      (
        score::RestraintNoeKnowledgeBased::GetDefaultScheme() == def_constr.GetScheme(),
        "Default def_constructor should return " + util::Format()( def_constr.GetScheme()) +
        " but instead returns " + util::Format()( score::RestraintNoeKnowledgeBased::GetDefaultScheme())
      );

      // test Clone constructor
      BCL_MessageStd( "test Clone def_constructor");
      const util::ShPtr< score::RestraintNoeKnowledgeBased> clone_constr( def_constr.Clone());
      BCL_Example_Check
      (
        def_constr.GetScheme() == clone_constr->GetScheme(),
        "Clone def_constructor should return " + util::Format()( def_constr.GetScheme()) +
        " but instead returned " + util::Format()( clone_constr->GetScheme())
      );

    /////////////////
    // data access //
    /////////////////

      // check GetClassIdentifier
      BCL_MessageStd( "test GetClassIdentifier function");
      const std::string correct_static_class_name( "bcl::score::RestraintNoeKnowledgeBased");
      BCL_Example_Check
      (
        GetStaticClassName< score::RestraintNoeKnowledgeBased>() == clone_constr->GetClassIdentifier() &&
        GetStaticClassName< score::RestraintNoeKnowledgeBased>() == correct_static_class_name,
        "GetClassIdentifier gives " + clone_constr->GetClassIdentifier() + " but should give " +
        correct_static_class_name
      );

      // test GetScheme
      const std::string correct_scheme( "noe_kb");
      BCL_ExampleCheck( def_constr.GetScheme(), correct_scheme);

    ///////////////
    // operators //
    ///////////////

      // test () operator for good restraint
      BCL_ExampleIndirectCheckWithinTolerance
      (
        def_constr( restraint_a.GenerateAssignment( protein_model)),
        -0.849824,
        0.0001,
        "operator () with good restraint"
      );

      // test () operator for ok restraint
      BCL_ExampleIndirectCheckWithinAbsTolerance
      (
        def_constr( restraint_b.GenerateAssignment( protein_model)),
        0.0,
        1e-10,
        "operator() with ok restraint"
      );

      // test () operator for bad restraint
      BCL_ExampleIndirectCheckWithinAbsTolerance
      (
        def_constr( restraint_c.GenerateAssignment( protein_model)),
        0.0,
        1e-10,
        "operator() with bad restraint"
      );

      // if message level is debug
      if( util::GetMessenger().IsSmallerEqualCurrentMessageLevel( util::Message::e_Debug))
      {
        // get the data
        const storage::Map
        <
          storage::Pair< biol::AtomType, size_t>,
          math::CubicSplineDamped
        > histogram_data( def_constr.GetNOEHistogram());

        // initialize splines
        util::SiPtrVector< const math::FunctionInterfaceSerializable< double, double> > sc_splines;
        storage::Vector< std::string> sc_descriptors;
        double sc_min( 1000);
        double sc_max( -1000);
        size_t sc_bin( 0);
        util::SiPtrVector< const math::FunctionInterfaceSerializable< double, double> > h_splines;
        storage::Vector< std::string> h_descriptors;
        double h_min( 1000);
        double h_max( -1000);
        size_t h_bin( 0);
        util::SiPtrVector< const math::FunctionInterfaceSerializable< double, double> > ha_splines;
        storage::Vector< std::string> ha_descriptors;
        double ha_min( 1000);
        double ha_max( -1000);
        size_t ha_bin( 0);

        // iterate over the data
        for
        (
          auto map_itr( histogram_data.Begin()), map_itr_end( histogram_data.End());
          map_itr != map_itr_end;
          ++map_itr
        )
        {
          // get the atom type
          const biol::AtomType &atom_type( map_itr->first.First());

          // store the spline depending on the atom type
          if( atom_type == biol::GetAtomTypes().CB)
          {
            sc_splines.PushBack( util::ToSiPtr( map_itr->second));
            sc_descriptors.PushBack( util::Format()( map_itr->first.Second()));
            sc_bin = map_itr->second.GetXValues().GetSize();
            sc_min = std::min( sc_min, map_itr->second.GetXValues().First());
            sc_max = std::max( sc_max, map_itr->second.GetXValues().Last());
          }
          else if( atom_type == biol::GetAtomTypes().H)
          {
            h_splines.PushBack( util::ToSiPtr( map_itr->second));
            h_descriptors.PushBack( util::Format()( map_itr->first.Second()));
            h_bin = map_itr->second.GetXValues().GetSize();
            h_min = std::min( h_min, map_itr->second.GetXValues().First());
            h_max = std::max( h_max, map_itr->second.GetXValues().Last());
          }
          else if( atom_type == biol::GetAtomTypes().HA)
          {
            ha_splines.PushBack( util::ToSiPtr( map_itr->second));
            ha_descriptors.PushBack( util::Format()( map_itr->first.Second()));
            ha_bin = map_itr->second.GetXValues().GetSize();
            ha_min = std::min( ha_min, map_itr->second.GetXValues().First());
            ha_max = std::max( ha_max, map_itr->second.GetXValues().Last());
          }
        }

        // write sc heatmap
        io::OFStream write;
        BCL_ExampleMustOpenOutputFile( write, "noe_kb_sc.gnuplot");
        math::GnuplotHeatmap sc_heatmap;
        sc_heatmap.SetFromFunctions
        (
          util::SiPtrVector< const math::FunctionInterfaceSerializable< double, double> >( sc_splines),
          sc_bin,
          sc_min,
          ( sc_max - sc_min) / double( sc_bin),
          false,
          false,
          false
        );
        sc_heatmap.SetTicsY( sc_descriptors, true, 1);
        sc_heatmap.SetRotationXTics( 90);
        sc_heatmap.SetPixelAndRatio( 2000, 800, -1);
        sc_heatmap.SetTitleAndLabel
        (
          "NOE knowledge based potential, side chain <--> side chain",
          "H-H distance - CB-CB distance [" + math::GnuplotHeatmap::s_AngstromSymbolGnuplot + "]",
          "# of bonds from CB",
          "Energy (BCLEUs)"
        );
        sc_heatmap.SetFont( "arialbd", 2);
        sc_heatmap.SetMinMaxZ( -1.0, 0.0);
        sc_heatmap.SetFilename( "noe_kb_sc");
        sc_heatmap.WriteScript( write);
        io::File::CloseClearFStream( write);

        // write h heatmap
        BCL_ExampleMustOpenOutputFile( write, "noe_kb_h.gnuplot");
        math::GnuplotHeatmap h_heatmap;
        h_heatmap.SetFromFunctions
        (
          util::SiPtrVector< const math::FunctionInterfaceSerializable< double, double> >( h_splines),
          ( h_max - h_min) / h_bin,
          h_min,
          h_bin,
          false,
          false,
          false
        );
        h_heatmap.SetTicsY( h_descriptors, true, 1);
        h_heatmap.SetRotationXTics( 90);
        h_heatmap.SetPixelAndRatio( 2000, 400, -1);
        h_heatmap.SetTitleAndLabel
        (
          "NOE knowledge based potential, H <--> side chain",
          "H-H distance - CB-CB distance [" + math::GnuplotHeatmap::s_AngstromSymbolGnuplot + "]",
          "# of bonds from CB",
          "Energy (BCLEUs)"
        );
        h_heatmap.SetFont( "arialbd", 2);
        h_heatmap.SetMinMaxZ( -1.0, 0.0);
        h_heatmap.SetFilename( "noe_kb_h");
        h_heatmap.WriteScript( write);
        io::File::CloseClearFStream( write);

        // write ha heatmap
        BCL_ExampleMustOpenOutputFile( write, "noe_kb_ha.gnuplot");
        math::GnuplotHeatmap ha_heatmap;
        ha_heatmap.SetFromFunctions
        (
          util::SiPtrVector< const math::FunctionInterfaceSerializable< double, double> >( ha_splines),
          ( ha_max - ha_min) / ha_bin,
          ha_min,
          ha_bin,
          false,
          false,
          false
        );
        ha_heatmap.SetTicsY( ha_descriptors, true, 1);
        ha_heatmap.SetRotationXTics( 90);
        ha_heatmap.SetPixelAndRatio( 2000, 400, -1);
        ha_heatmap.SetTitleAndLabel
        (
          "NOE knowledge based potential, HA <--> side chain",
          "H-H distance - CB-CB distance [" + math::GnuplotHeatmap::s_AngstromSymbolGnuplot + "]",
          "# of bonds from CB",
          "Energy (BCLEUs)"
        );
        ha_heatmap.SetFont( "arialbd", 2);
        ha_heatmap.SetMinMaxZ( -1.0, 0.0);
        ha_heatmap.SetFilename( "noe_kb_ha");
        ha_heatmap.WriteScript( write);
        io::File::CloseClearFStream( write);
      }

    //////////////////////
    // input and output //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleScoreRestraintNoeKnowledgeBased

  const ExampleClass::EnumType ExampleScoreRestraintNoeKnowledgeBased::s_Instance
  (
    GetExamples().AddEnum( ExampleScoreRestraintNoeKnowledgeBased())
  );

} // namespace bcl
