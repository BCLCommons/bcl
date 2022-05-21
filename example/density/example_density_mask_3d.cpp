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
#include "density/bcl_density_mask_3d.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "biol/bcl_biol_atom.h"
#include "density/bcl_density_map.h"
#include "density/bcl_density_simulators.h"
#include "math/bcl_math_running_average_sd.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_density_mask_3d.cpp
  //!
  //! @author bitterd
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //! @remarks need to be reviewed
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleDensityMask3d :
    public ExampleInterface
  {
  public:

    ExampleDensityMask3d *Clone() const
    { return new ExampleDensityMask3d( *this);}

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
      const std::string example_pdb( AddExampleInputPathToFilename( e_Biology, "1ubi.pdb"));

    /////////////////
    // preparation //
    /////////////////

      // create pdb model
      storage::Map< biol::SSType, size_t> sse_min_size;
      sse_min_size[ biol::GetSSTypes().HELIX] = 0;
      sse_min_size[ biol::GetSSTypes().STRAND] = 0;
      sse_min_size[ biol::GetSSTypes().COIL] = 0;
      assemble::ProteinModel protein( Proteins::GetModel( example_pdb, biol::GetAAClasses().e_AAComplete, sse_min_size));

      // set parameters
      const double resolution( 6.6), voxelsize( 2.2);
      const util::ShPtr< density::SimulateInterface> simulator
      (
        density::GetSimulators().CreateSimulator
        (
          density::GetSimulators().e_Gaussian, linal::Vector3D( voxelsize), resolution
        )
      );
      const double masking_distance( 5.0);

      // density map calculated from SimpleAtoms, with target resolution, voxelsize( usually sigma = resolution/3) and smoothingkernel of choice
      const density::Map density_sim( simulator->operator ()( protein.GetAtoms()));
      density::Map density_exp( density_sim);
      density_exp.AddNoise( random::GetGlobalRandom(), double( 0.0), 1.5 * density_sim.GetRmsd());

      // create mask for every three consecutive residues
      const util::SiPtrVector< const biol::AABase> residues( protein.GetAminoAcids());
      const size_t nr_residues( residues.GetSize());

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // create map of pointers to atoms and corresponding masks
      storage::Map< util::SiPtr< const biol::AABase>, density::Mask3d> density_masks;
      for( size_t res( 1); res < nr_residues - 1; ++res)
      {
        // get set of 3 amino acids centered around aa of interest
        util::SiPtrVector< const linal::Vector3D> coords;
        coords.Append( residues( res - 1)->GetAtomCoordinates());
        coords.Append( residues( res    )->GetAtomCoordinates());
        coords.Append( residues( res + 1)->GetAtomCoordinates());

        // insert constructed mask into map
        density_masks[ residues( res)] =
          density::Mask3d( coords, masking_distance, density_exp.GetCellWidth(), density_exp.GetOrigin());
      }

      BCL_Example_Check
      (
        density_masks.GetSize() == nr_residues - 2,
        "map should contain an entry for every residue but the first and the last one ( " +
          util::Format()( nr_residues - 2) + " amino acids), but contains " +
          util::Format()( density_masks.GetSize()) + "entries"
      );

    /////////////////
    // data access //
    /////////////////

      // define bool to keep track of success
      bool success( true);

      // iterate over residues and their masks
      for
      (
        storage::Map< util::SiPtr< const biol::AABase>, density::Mask3d>::const_iterator
          itr( density_masks.Begin()), itr_end( density_masks.End());
        itr != itr_end;
        ++itr
      )
      {
        const math::RunningAverageSD< double> statistics_sim( itr->second.CalculateMeanSD( density_sim));
        const math::RunningAverageSD< double> statistics_exp( itr->second.CalculateMeanSD( density_exp));
        const double correlation( itr->second.CrossCorrelationCoefficient( density_exp, density_sim));

        BCL_MessageStd
        (
          "mean for mask res: " + itr->first->GetIdentification() +
            " statistics_exp: " + util::Format()( statistics_exp.GetAverage())
            + " statistics_sim: " + util::Format()( statistics_sim.GetAverage())
        );
        BCL_MessageStd
        (
          "sd for mask res: " + itr->first->GetIdentification() +
            " statistics_exp: " + util::Format()( statistics_exp.GetStandardDeviation()) +
            " statistics_sim: " + util::Format()( statistics_sim.GetStandardDeviation())
        );
        BCL_MessageStd
        (
          "CrossCorrelationCoefficient: " + itr->first->GetIdentification()
            + ":  " + util::Format()( correlation)
        );
        density::Mask3d mask_copy( itr->second);

        // test whether the copy was done correctly
        if
        (
          !(
            itr->second.GetIndex() == mask_copy.GetIndex()
            && itr->second.GetPosition() == mask_copy.GetPosition()
            && itr->second.GetGridSpacing() == mask_copy.GetGridSpacing()
            && itr->second.GetMask() == mask_copy.GetMask()
          )
        )
        {
          BCL_MessageStd
          (
            "Copying failed in residue" + itr->first->GetIdentification()
            + util::Format()( itr->second.GetIndex()) + "!=" + util::Format()( mask_copy.GetIndex()) + " or "
            + util::Format()( itr->second.GetPosition()) + "!=" + util::Format()( mask_copy.GetPosition()) + " or "
            + util::Format()( itr->second.GetGridSpacing()) + "!=" + util::Format()( mask_copy.GetGridSpacing()) + " or "
            + util::Format()( itr->second.GetMask()) + "!=" + util::Format()( mask_copy.GetMask())
          );
          success = false;
        }

        // do a clone
        util::ShPtr< density::Mask3d> mask_clone( itr->second.Clone());

        // test whether the clone was done correctly
        if
        (
          !(
            itr->second.GetIndex() == mask_clone->GetIndex()
            && itr->second.GetPosition() == mask_clone->GetPosition()
            && itr->second.GetGridSpacing() == mask_clone->GetGridSpacing()
            && itr->second.GetMask() == mask_clone->GetMask()
          )
        )
        {
          BCL_MessageStd
          (
            "Cloning failed in residue" + itr->first->GetIdentification()
            + util::Format()( itr->second.GetIndex()) + "!=" + util::Format()( mask_clone->GetIndex()) + " or "
            + util::Format()( itr->second.GetPosition()) + "!=" + util::Format()( mask_clone->GetPosition()) + " or "
            + util::Format()( itr->second.GetGridSpacing()) + "!=" + util::Format()( mask_clone->GetGridSpacing()) + " or "
            + util::Format()( itr->second.GetMask()) + "!=" + util::Format()( mask_clone->GetMask())
          );
          success = false;
        }

        BCL_MessageStd
        (
          "count of considered values for res: " + itr->first->GetIdentification()
            + " statistics_exp: " + util::Format()( statistics_exp.GetWeight())
            + " statistics_sim: " + util::Format()( statistics_sim.GetWeight())
        );
      }

      // test whether all checks were successful
      BCL_Example_Check( success, "Example has failed");

      // test GetStaticClassName and GetClassIdentifier function
      BCL_Example_Check
      (
        GetStaticClassName< density::Mask3d>() == "bcl::density::Mask3d"
          && density_masks.Begin()->second.GetClassIdentifier() == GetStaticClassName< density::Mask3d>(),
        "GetStaticClassName and GetClassIdentifier functions don't work correctly"
      );

    //////////////////////
    // input and output //
    //////////////////////

      // density mask to write
      const density::Mask3d &mask_to_write( density_masks[ residues( 1)]);
      WriteBCLObject( mask_to_write);

      // read mask
      density::Mask3d mask3d_read;
      ReadBCLObject( mask3d_read);

      BCL_Example_Check
      (
        mask_to_write.GetIndex()       == mask3d_read.GetIndex() &&
        mask_to_write.GetPosition()    == mask3d_read.GetPosition() &&
        mask_to_write.GetGridSpacing() == mask3d_read.GetGridSpacing(),
        "Read and write of mask failed for second residue" +
        util::Format()( mask_to_write.GetIndex())       + "!=" + util::Format()( mask3d_read.GetIndex()) + " or " +
        util::Format()( mask_to_write.GetPosition())    + "!=" + util::Format()( mask3d_read.GetPosition()) + " or " +
        util::Format()( mask_to_write.GetGridSpacing()) + "!=" + util::Format()( mask3d_read.GetGridSpacing()) + " or " +
        util::Format()( mask_to_write.GetMask())        + "!=" + util::Format()( mask3d_read.GetMask())
      );

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; // end ExampleDensityMask3d

  const ExampleClass::EnumType ExampleDensityMask3d::s_Instance
  (
    GetExamples().AddEnum( ExampleDensityMask3d())
  );

} // namespace bcl
