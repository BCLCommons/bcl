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
#include "density/bcl_density_map_cylindrical.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "biol/bcl_biol_atom.h"
#include "density/bcl_density_map.h"
#include "density/bcl_density_simulators.h"
#include "io/bcl_io_file.h"
#include "math/bcl_math_gnuplot_heatmap.h"
#include "math/bcl_math_histogram.h"
#include "math/bcl_math_histogram_2d.h"
#include "pdb/bcl_pdb_factory.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_density_map_cylindrical.cpp
  //!
  //! @author linders
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleDensityMapCylindrical :
    public ExampleInterface
  {
  public:

    ExampleDensityMapCylindrical *Clone() const
    { return new ExampleDensityMapCylindrical( *this);}

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

    /////////////////
    // preparation //
    /////////////////

      // initialize pdb filename
      //const std::string pdb_filename( AddExampleInputPathToFilename( e_Biology, "1IE9_first_helix.pdb"));
      const std::string pdb_filename( AddExampleInputPathToFilename( e_Biology, "1IE9.pdb"));
      //const std::string pdb_filename( AddExampleInputPathToFilename( e_Biology, "hex1_fit_rnd38_v3.pdb"));
      BCL_MessageStd( "reading pdb: " + pdb_filename);

      storage::Map< biol::SSType, size_t> complete_ssetype_min_size;
      complete_ssetype_min_size[ biol::GetSSTypes().HELIX] = 0;
      complete_ssetype_min_size[ biol::GetSSTypes().STRAND] = 0;
      complete_ssetype_min_size[ biol::GetSSTypes().COIL] = 0;
      const assemble::ProteinModel complete_model( Proteins::GetModel( pdb_filename, biol::GetAAClasses().e_AAComplete, complete_ssetype_min_size));
//      complete_model.SetToIdealConformation();
//      complete_model = *biol::AASideChainFactory( false, true).ProteinModelWithSideChains( complete_model);

      //read atoms from pdb
      io::IFStream read;

      //set parameters of that density map
      const double resolution( 6.9), voxelsize( 1.5);
      const util::ShPtr< density::SimulateInterface> simulator( density::GetSimulators().CreateSimulator( density::GetSimulators().e_Gaussian, linal::Vector3D( voxelsize), resolution));

      // density map calculated from SimpleAtoms, with target resolution, voxelsize( usually sigma = resolution/3)
      // and smoothingkernel of choice
      const density::Map density_map( simulator->operator ()( complete_model.GetAtoms()));

      // next steps are to extract a SSE body from the pdb
      // create factory
      pdb::Factory factory( biol::GetAAClasses().e_AAComplete);

      //build protein models from pdb sequence
      BCL_MessageStd( "building models from pdb chains and sse information");
      storage::Map< biol::SSType, size_t> ssetype_min_size;
      ssetype_min_size[ biol::GetSSTypes().HELIX] = 10;
      ssetype_min_size[ biol::GetSSTypes().STRAND] = 99;
      assemble::ProteinModel protein_model( Proteins::GetModel( pdb_filename, biol::GetAAClasses().e_AAComplete, ssetype_min_size));

      // initialize body from one of the SSEs in the protein model
      const assemble::SSEGeometry sse_body_first_helix( *protein_model.GetSSEs().FirstElement());
      BCL_MessageStd
      (
        "z-axis dimension of first SSE in 1IE9: " +
        util::Format()( sse_body_first_helix.GetExtent( coord::GetAxes().e_Z))
      );

      const util::SiPtrVector< const assemble::SSEGeometryInterface> temp1( protein_model.GetSSEs());

      // initialize SiPtrList of all the bodies of helices in the protein
      const util::SiPtrList< const assemble::SSEGeometryInterface> list_of_sses
      (
        temp1.Begin(),
        temp1.End()
      );

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // default constructor
      const density::MapCylindrical density_map_cylindrical_default;

      const double height_resolution( 1.5);
      const double radius_resolution( 1.0);
      const size_t number_wedges( 36);
      const double lower_radius ( 2);
      const double upper_radius( 6);

      // constructor from (Cartesian) density map and body
      const density::MapCylindrical density_map_cylindrical_first_helix
      (
        sse_body_first_helix,
        density_map,
        height_resolution,
        radius_resolution,
        number_wedges,
        upper_radius
      );

      // test constructor from list of bodies and spline indirectly by calling static function
      // MapCylindrical::CalculateCylindricalMaps
      const storage::List< density::MapCylindrical> list_of_cylindrical_density_maps
      (
        density::MapCylindrical::CalculateCylindricalMaps
        (
          list_of_sses,
          density_map,
          height_resolution,
          radius_resolution,
          number_wedges,
          upper_radius
        )
      );

    /////////////////
    // data access //
    /////////////////

      BCL_MessageStd
      (
        "dimensions of cylindrical density map: " +
        util::Format()( density_map_cylindrical_first_helix.GetDimensions())
      );

      // test operator
      const double expected_intensity_at_center( 8.85472);
      const double expected_intensity_at_edge( 5.88442);
      const double intensity_at_center( density_map_cylindrical_first_helix( 9, 0, 0));
      const double intensity_at_edge
      (
        density_map_cylindrical_first_helix( 9, size_t( upper_radius / radius_resolution) - 1, 0)
      );
      BCL_MessageStd
      (
        "intensity at center, intensity at edge: " +
        util::Format()( intensity_at_center) + ", " + util::Format()( intensity_at_edge)
      );
      BCL_ExampleCheckWithinTolerance( intensity_at_center, expected_intensity_at_center, 0.0001);
      BCL_ExampleCheckWithinTolerance( intensity_at_edge, expected_intensity_at_edge, 0.0001);

    ////////////////
    // operations //
    ////////////////

      size_t helix_counter( 1);
      // iterate over list of cylindrical density maps and output 2D profiles / gnuplot heatmap files
      for
      (
        storage::List< density::MapCylindrical>::const_iterator itr_density_map
          ( list_of_cylindrical_density_maps.Begin()),
          itr_density_map_end( list_of_cylindrical_density_maps.End());
        itr_density_map != itr_density_map_end;
        ++itr_density_map
      )
      {

        // calculate TwoDProfileHeightAngle
        const math::Histogram2D height_angle_profile( itr_density_map->TwoDProfileHeightAngle( lower_radius, upper_radius));

        // calculate TwoDProfileHeightRadius
        const math::Histogram2D height_radius_profile( itr_density_map->TwoDProfileHeightRadius( lower_radius, upper_radius));

        // calculate TwoDProfileHeightRadius
        const math::Histogram2D radius_angle_profile( itr_density_map->TwoDProfileRadiusAngle( lower_radius, upper_radius));

        //calculate OneDProfileHeight
        const math::Histogram height_profile( itr_density_map->OneDProfileHeight( lower_radius, upper_radius));

        // calculate OneDProfileRadius
        const math::Histogram radius_profile( itr_density_map->OneDProfileRadius( lower_radius, upper_radius));

        //calculate OneDProfileHeight
        const math::Histogram angle_profile( itr_density_map->OneDProfileAngle( lower_radius, upper_radius));

        if( util::GetMessenger().IsSmallerEqualCurrentMessageLevel( util::Message::e_Debug))
        {
          // output the gnuplot heatmap files
          io::OFStream write;

          BCL_ExampleMustOpenOutputFile( write, util::Format()( helix_counter) + "_two_d_profile_height_radius.gnuplot");
          math::GnuplotHeatmap heatmap_height_radius_profile;
          heatmap_height_radius_profile.SetFromHistogram( height_radius_profile, false, false);
          heatmap_height_radius_profile.SetTitleAndLabel( util::Format()( helix_counter) + "_2D_profile_height_radius", "height in rod [" + math::GnuplotHeatmap::s_AngstromSymbolGnuplot + "]", "distance from center [" + math::GnuplotHeatmap::s_AngstromSymbolGnuplot + "]", "density sum (angle)");
          heatmap_height_radius_profile.SetFont( "arialbd", 16);
          heatmap_height_radius_profile.SetRotationXTics( 90.0);
          heatmap_height_radius_profile.SetFilename( util::Format()( helix_counter) + "_two_d_profile_height_radius");
          heatmap_height_radius_profile.WriteScript( write);
          io::File::CloseClearFStream( write);

          BCL_ExampleMustOpenOutputFile( write, util::Format()( helix_counter) + "_two_d_profile_height_angle.gnuplot");
          math::GnuplotHeatmap heatmap_height_angle_profile;
          heatmap_height_angle_profile.SetFromHistogram( height_angle_profile, false, false);
          heatmap_height_angle_profile.SetTitleAndLabel( util::Format()( helix_counter) + "_2D_profile_height_angle", "height in rod [" + math::GnuplotHeatmap::s_AngstromSymbolGnuplot + "]", "angle [rad]", "density sum (height)");
          heatmap_height_angle_profile.SetFont( "arialbd", 16);
          heatmap_height_angle_profile.SetRotationXTics( 90.0);
          heatmap_height_angle_profile.SetFilename( util::Format()( helix_counter) + "_two_d_profile_height_angle");
          heatmap_height_angle_profile.WriteScript( write);
          io::File::CloseClearFStream( write);

          BCL_ExampleMustOpenOutputFile( write, util::Format()( helix_counter) + "_two_d_profile_radius_angle.gnuplot");
          math::GnuplotHeatmap heatmap_radius_angle_profile;
          heatmap_radius_angle_profile.SetFromHistogram( radius_angle_profile, false, false);
          heatmap_radius_angle_profile.SetTitleAndLabel( util::Format()( helix_counter) + "_2D_profile_radius_angle", "distance from center [" + math::GnuplotHeatmap::s_AngstromSymbolGnuplot + "]", "angle [rad]", "density sum (height)");
          heatmap_radius_angle_profile.SetFont( "arialbd", 16);
          heatmap_radius_angle_profile.SetFilename( util::Format()( helix_counter) + "_two_d_profile_radius_angle");
          heatmap_radius_angle_profile.WriteScript( write);
          io::File::CloseClearFStream( write);

          BCL_ExampleMustOpenOutputFile( write, util::Format()( helix_counter) + "_one_d_profile_height.gnuplot");
          math::GnuplotHeatmap heatmap_height_profile;
          heatmap_height_profile.SetFromHistogram( height_profile, false, false);
          heatmap_height_profile.SetTitleAndLabel( util::Format()( helix_counter) + "_1D_profile_height", "height in rod [" + math::GnuplotHeatmap::s_AngstromSymbolGnuplot + "]", "", "density sum (angle and radius)");
          heatmap_height_profile.SetFont( "arialbd", 16);
          heatmap_height_profile.SetRotationXTics( 90.0);
          heatmap_height_profile.SetFilename( util::Format()( helix_counter) + "_one_d_profile_height");
          heatmap_height_profile.WriteScript( write);
          io::File::CloseClearFStream( write);

          BCL_ExampleMustOpenOutputFile( write, util::Format()( helix_counter) + "_one_d_profile_radius.gnuplot");
          math::GnuplotHeatmap heatmap_radius_profile;
          heatmap_radius_profile.SetFromHistogram( radius_profile, false, false);
          heatmap_radius_profile.SetTitleAndLabel( util::Format()( helix_counter) + "_1D_profile_radius", "distance from center [" + math::GnuplotHeatmap::s_AngstromSymbolGnuplot + "]", "", "density sum (angle and height)");
          heatmap_radius_profile.SetFont( "arialbd", 16);
          heatmap_radius_profile.SetFilename( util::Format()( helix_counter) + "_one_d_profile_radius");
          heatmap_radius_profile.WriteScript( write);
          io::File::CloseClearFStream( write);

          BCL_ExampleMustOpenOutputFile( write, util::Format()( helix_counter) + "_one_d_profile_angle.gnuplot");
          math::GnuplotHeatmap heatmap_angle_profile;
          heatmap_angle_profile.SetFromHistogram( angle_profile, false, false);
          heatmap_angle_profile.SetTitleAndLabel( util::Format()( helix_counter) + "_1D_profile_angle", "angle [rad]", "", "density sum (radius and height)");
          heatmap_angle_profile.SetFont( "arialbd", 16);
          heatmap_angle_profile.SetRotationXTics( 90.0);
          heatmap_angle_profile.SetFilename( util::Format()( helix_counter) + "_one_d_profile_angle");
          heatmap_angle_profile.WriteScript( write);
          io::File::CloseClearFStream( write);

          ++helix_counter;
        }
      }

    //////////////////////
    // input and output //
    //////////////////////

      BCL_MessageStd( "testing read and write functionalities for density::MapCylindrical");
      // open output file
      BCL_MessageVrb( "write object");
      WriteBCLObject( density_map_cylindrical_first_helix);
      BCL_MessageVrb( "read object");
      density::MapCylindrical density_map_read;
      ReadBCLObject( density_map_read);

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; // end ExampleDensityMapCylindrical

  const ExampleClass::EnumType ExampleDensityMapCylindrical::s_Instance
  (
    GetExamples().AddEnum( ExampleDensityMapCylindrical())
  );

} // namespace bcl
