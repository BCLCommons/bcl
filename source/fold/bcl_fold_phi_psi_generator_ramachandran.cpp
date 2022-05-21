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
#include "fold/bcl_fold_phi_psi_generator_ramachandran.h"

// includes from bcl - sorted alphabetically
#include "fold/bcl_fold_mutation_residue.h"
#include "io/bcl_io_serialization.h"
#include "util/bcl_util_enumerated.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {
    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> PhiPsiGeneratorRamachandran::s_Instance
    (
      util::Enumerated< math::FunctionInterfaceSerializable< MutationResidue, storage::VectorND< 2, double> > >::AddInstance
      (
        new PhiPsiGeneratorRamachandran()
      )
    );

    //! @brief get default instance
    //! @return default instance behind a pointer
    const math::FunctionInterfaceSerializable< MutationResidue, storage::VectorND< 2, double> > &
    PhiPsiGeneratorRamachandran::GetDefaultInstance()
    {
      // initialize static instance
      static const PhiPsiGeneratorRamachandran s_instance( biol::Ramachandran::GetDefaultInstance());

      // end
      return s_instance;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    PhiPsiGeneratorRamachandran::PhiPsiGeneratorRamachandran() :
      m_Ramachandran( biol::Ramachandran())
    {
    }

    //! @brief constructor from a ramachandran object
    //! @param RAMACHANDRAN instance of Ramachandran to be used
    PhiPsiGeneratorRamachandran::PhiPsiGeneratorRamachandran( const biol::Ramachandran &RAMACHANDRAN) :
       m_Ramachandran( RAMACHANDRAN)
    {
    }

    //! @brief Clone function
    //! @return pointer to new PhiPsiGeneratorRamachandran
    PhiPsiGeneratorRamachandran *PhiPsiGeneratorRamachandran::Clone() const
    {
      return new PhiPsiGeneratorRamachandran( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &PhiPsiGeneratorRamachandran::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get the name of the object when used in a dynamic context
    //! @return the name of the object when used in a dynamic context
    const std::string &PhiPsiGeneratorRamachandran::GetAlias() const
    {
      static const std::string s_alias( "PhiPsiGeneratorRamachandran");
      return s_alias;
    }

    //! @brief returns parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer PhiPsiGeneratorRamachandran::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "Generate random pairs of phi psi angles.");
      serializer.AddInitializer
      (
        "distribution",
        "phi psi distribution from which to generate the angles",
        io::Serialization::GetAgent( &m_Ramachandran),
        "Ramachandran"
      );

      return serializer;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief operator taking MutationResidue and returning storage::VectorND< 2, double> with phi psi, respectively
    //! @param RESIDUE MutationResidue whose phi or psi will be changed
    //! @return phi or psi value, respectively, randomly selected based on a ramachandran distribution
    storage::VectorND< 2, double> PhiPsiGeneratorRamachandran::operator()( const MutationResidue &RESIDUE) const
    {
      BCL_Assert( RESIDUE.GetMutationResidue().IsDefined(), "mutation residue is not defined");
      // get random phi and psi
      const storage::VectorND< 2, double> random_phi_psi
      (
        m_Ramachandran->GetRandomPhiPsi( RESIDUE.GetMutationResidue()->GetType())
      );

      // return "random_phi_psi"
      return random_phi_psi;
    }

  //////////////////////
  // input and output //
  //////////////////////

  } // namespace fold
} // namespace bcl
