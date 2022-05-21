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
#include "chemistry/bcl_chemistry_mutate_molecule_generic.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_atom_clash_score.h"
#include "chemistry/bcl_chemistry_configurational_bond_type_data.h"
#include "chemistry/bcl_chemistry_conformation_comparison_property_field_correlation.h"
#include "chemistry/bcl_chemistry_conformation_interface.h"
#include "chemistry/bcl_chemistry_mutate_fragment.h"
#include "descriptor/bcl_descriptor_constants.h"
#include "graph/bcl_graph_undirected_edge.h"
#include "linal/bcl_linal_vector_3d_operations.h"
#include "mc/bcl_mc_temperature_accepted.h"
#include "quality/bcl_quality_rmsd.h"
#include "storage/bcl_storage_pair.h"
// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

    // add the interface to the set of known implementations
    const util::SiPtr< const util::ObjectInterface> MutateMoleculeGeneric::s_Instance
    (
      GetObjectInstances().AddInstance( new MutateMoleculeGeneric())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    MutateMoleculeGeneric::MutateMoleculeGeneric()
    {
    }

    //! @brief movie constructor
    MutateMoleculeGeneric::MutateMoleculeGeneric
    (
      const std::string &MOVIE_FILENAME
    ) :
        m_Movie( MOVIE_FILENAME)
    {
    }

    //! virtual copy constructor
    MutateMoleculeGeneric *MutateMoleculeGeneric::Clone() const
    {
      return new MutateMoleculeGeneric( *this);
    }

  /////////////////
  // data access // 
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &MutateMoleculeGeneric::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns the name used for this class in an object data label
    //! @return the name used for this class in an object data label
    const std::string &MutateMoleculeGeneric::GetAlias() const
    {
      static std::string s_name( "MoleculeMutateGeneric");
      return s_name;
    }

  /////////////////
  //  operations //
  /////////////////

    //! @brief align two small molecule objects and find the property RMSD
    //! @param MOLECULE - molecule to be fit into pocket
    //! @param POCKET - pocket into which MOLECULE is being geometrically fit
    //! @return the molecule in its new pose relative to the static pocket
    math::MutateResult< FragmentComplete> MutateMoleculeGeneric::operator()
    (
      const FragmentComplete &MOLECULE
    ) const
    {
        FragmentComplete mol( MOLECULE);
        double prob( random::GetGlobalRandom().Random< double>( 1));

        // rotate molecule randomly with magnitude between 0 and 5 degrees
        if( prob < 0.25)
          {
            math::RotationMatrix3D rot;
            rot.SetRand( 0.08);
            linal::Vector3D mol_centered( mol.GetCenter());
            mol.Translate( -mol_centered);
            mol.Rotate( rot);
            mol.Translate( mol_centered);
          }

        // translate molecule a random distance between 0 and 1 angstroms
        else if( prob < 0.50)
          {
            const double mag( random::GetGlobalRandom().Random( 1.0));
            linal::Vector3D trans( mag, 0.0, 0.0);
            mol.Translate( trans.Rotate( math::RotationMatrix3D().SetRand()));
          }

        // rotate molecule randomly with magnitude between 0 and 180 degrees
        else if( prob < 0.75)
          {
            math::RotationMatrix3D rot;
            rot.SetRand();
            linal::Vector3D mol_centered( mol.GetCenter());
            mol.Translate( -mol_centered);
            mol.Rotate( rot);
            mol.Translate( mol_centered);
          }
        else
          {
            const storage::Vector< float> flip_angles( storage::Vector< float>::Create( 90.0, 180.0));
            linal::Vector3D axis;
            axis.SetRandomTranslation( 1.0);
            math::RotationMatrix3D flip( axis, flip_angles( size_t( ( prob - 0.75) / 0.125)));
            linal::Vector3D mol_centered( mol.GetCenter());
            mol.Translate( -mol_centered);
            mol.Rotate( flip);
            mol.Translate( mol_centered);
          }

        //And now output the transformed molecule with the rest of the mutate info
        util::ShPtr< FragmentComplete> mutant( new FragmentComplete( mol));
        mutant->GetCacheMap() = MOLECULE.GetCacheMap();
        math::MutateResult< FragmentComplete> mutate_result( mutant, *this);

        // Outputs each transformation to an .sdf file
        // Primarily used for debug purposes
        if( !m_Movie.empty())
          {
          io::OFStream out;
          io::File::MustOpenOFStream( out, m_Movie + ".sdf", std::ios::app);
          mutant->WriteMDL( out);
          }

        return mutate_result;

    }

     //! @brief return parameters for member data that are set up from the labels
     //! @return parameters for member data that are set up from the labels
    io::Serializer MutateMoleculeGeneric::GetSerializer() const
    {
      io::Serializer member_data;
      member_data.SetClassDescription( "Orients a small molecule in a pocket cavity by minimizing geometric overlap of matched atoms");
      return member_data;

      member_data.AddOptionalInitializer
      (
        "output movie",
        "output each step of the attempted orientation of the ligand into the pocket - WARNING - large file size easily accumulates",
        io::Serialization::GetAgent( &m_Movie)
      );
    }
    //! @brief Set the members of this property from the given LABEL
    //! @param LABEL the label to parse
    //! @param ERR_STREAM stream to write out errors to
    bool MutateMoleculeGeneric::ReadInitializerSuccessHook
    (
      const util::ObjectDataLabel &LABEL,
      std::ostream &ERR_STREAM
    )
    {
      return math::MutateInterface< FragmentComplete>::ReadInitializerSuccessHook( LABEL, ERR_STREAM);
    }

  } // namespace chemistry
} // namespace bcl

