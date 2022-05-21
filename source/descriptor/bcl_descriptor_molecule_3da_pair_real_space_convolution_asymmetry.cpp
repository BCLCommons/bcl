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
#include "descriptor/bcl_descriptor_molecule_3da_pair_real_space_convolution_asymmetry.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_aa_side_chain_factory.h"
#include "chemistry/bcl_chemistry_fragment_ensemble.h"
#include "command/bcl_command_command_state.h"
#include "descriptor/bcl_descriptor_molecule_3da_smooth.h"
#include "io/bcl_io_ifstream.h"
#include "iterate/bcl_iterate_reflecting.h"
#include "math/bcl_math_running_min_max.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {

    // add each of the possible instances to the enumerated instances
    const util::SiPtr< const util::ObjectInterface> Molecule3DAPairRealSpaceConvolutionAsymmetry::s_Instance
    (
      util::Enumerated< Base< chemistry::AtomConformationalInterface, float> >::AddInstance
      (
        new Molecule3DAPairRealSpaceConvolutionAsymmetry()
      )
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    Molecule3DAPairRealSpaceConvolutionAsymmetry::Molecule3DAPairRealSpaceConvolutionAsymmetry() :
        Molecule3DAPairConvolutionAsymmetry()
    {
      if( !command::CommandState::IsInStaticInitialization())
      {

        BCL_Assert
        (
          ReadInitializerSuccessHook( util::ObjectDataLabel(), util::GetLogger()),
          "Failed to create " + GetClassIdentifier()
        );
      }
    }

    //! @brief constructor from number of steps, and mapped atom property
    Molecule3DAPairRealSpaceConvolutionAsymmetry::Molecule3DAPairRealSpaceConvolutionAsymmetry
    (
      const CheminfoProperty &ATOM_PROPERTY_A,
      const CheminfoProperty &ATOM_PROPERTY_B,
      const size_t NUMBER_STEPS,
      const float STEP_SIZE,
      const size_t WINDOW_SIZE,
      const float CUTOFF
    ) :
      Molecule3DAPairConvolutionAsymmetry( ATOM_PROPERTY_A, ATOM_PROPERTY_B, NUMBER_STEPS, STEP_SIZE, WINDOW_SIZE),
      m_Cutoff( CUTOFF),
      m_VoxelGridPocket( m_Cutoff)
    {
      BCL_Assert
      (
        ReadInitializerSuccessHook( util::ObjectDataLabel(), util::GetLogger()),
        "Failed to create " + GetClassIdentifier()
      );
    }

    //! @brief virtual copy constructor
    //! @return pointer to new Molecule3DAPairRealSpaceConvolutionAsymmetry
    Molecule3DAPairRealSpaceConvolutionAsymmetry *Molecule3DAPairRealSpaceConvolutionAsymmetry::Clone() const
    {
      return new Molecule3DAPairRealSpaceConvolutionAsymmetry( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief get name of the current class
    //! @return name of the class
    const std::string &Molecule3DAPairRealSpaceConvolutionAsymmetry::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the name of the property without any parameters
    //! @return name of the property as string
    const std::string &Molecule3DAPairRealSpaceConvolutionAsymmetry::GetAlias() const
    {
      static const std::string s_name( "3DAPairRealSpaceConvolutionAsymmetry");
      return s_name;
    }

  ////////////////
  // operations //
  ////////////////

    void Molecule3DAPairRealSpaceConvolutionAsymmetry::ReduceProteinPocket( chemistry::FragmentComplete &POCKET)
    {
      //Setup voxel grid for protein
      m_VoxelGridPocket.SetObjects( util::SiPtrVector< const chemistry::AtomConformationalInterface>( POCKET.GetAtomsIterator(), POCKET.GetAtomsIterator().End()));
      storage::Set< size_t> atom_indices;
      for
      (
          auto itr_m( GetCurrentObject()->GetIterator());
          itr_m.NotAtEnd();
          ++itr_m
      )
      {
        const chemistry::AtomConformationalInterface &atom_m( *itr_m);

        // Find neighbors_p of atom_m in pocket
        auto neighbors_p( m_VoxelGridPocket.GetNeighbors( atom_m.GetPosition(), m_Cutoff));

        // Get atominfo from every neighbor so we can make a new FragmentComplete
        for( auto itr_neighbors( neighbors_p.Begin()), itr_neighbors_end( neighbors_p.End()); itr_neighbors != itr_neighbors_end; ++itr_neighbors)
        {
          atom_indices.Insert( POCKET.GetAtomIndex( *itr_neighbors->First()));
        }
      }
      chemistry::AtomVector< chemistry::AtomComplete> atoms( POCKET.GetAtomVector());
      storage::Vector< size_t> atom_vector( atom_indices.Begin(), atom_indices.End());
      atoms.Reorder( atom_vector);
      POCKET = chemistry::FragmentComplete( atoms, POCKET.GetName());
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer Molecule3DAPairRealSpaceConvolutionAsymmetry::GetSerializer() const
    {
      io::Serializer parameters( Molecule3DAPairConvolutionAsymmetry::GetSerializer());
      parameters.SetClassDescription
      (
        "Relates the autocorrelations of two independent molecules"
      );

      // set the cutoff for identifying neighbor atoms
      parameters.AddInitializer
      (
        "cutoff",
        "distance cutoff (in Angstroms) for neighbor atom determination",
        io::Serialization::GetAgent( &m_Cutoff),
        "7.0"
      );
      return parameters;
    }

    //! @brief Set the members of this property from the given LABEL
    //! @param LABEL the label to parse
    //! @param ERR_STREAM stream to write out errors to
    bool Molecule3DAPairRealSpaceConvolutionAsymmetry::ReadInitializerSuccessHook
    (
      const util::ObjectDataLabel &LABEL,
      std::ostream &ERR_STREAM
    )
    {
      bool base( Molecule3DAPairConvolutionAsymmetry::ReadInitializerSuccessHook( LABEL, ERR_STREAM));
      m_VoxelGridPocket = chemistry::VoxelGridAtom( m_Cutoff);
      return base;
    }

    //! @brief hook that derived classes can override to add behavior after every time SetObject is called
    void Molecule3DAPairRealSpaceConvolutionAsymmetry::SetObjectHook()
    {
      // Reduce the protein pocket and then perform calculate in parent class
      Molecule3DAPairConvolutionAsymmetry::SetObjectHook();
    }

  } // namespace descriptor
} // namespace bcl
