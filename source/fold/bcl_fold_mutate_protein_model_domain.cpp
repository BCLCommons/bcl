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
#include "fold/bcl_fold_mutate_protein_model_domain.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_domain.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "math/bcl_math_mutate_move_wrapper.h"
#include "math/bcl_math_mutate_result.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> MutateProteinModelDomain::s_Instance
    (
      GetObjectInstances().AddInstance( new MutateProteinModelDomain())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    MutateProteinModelDomain::MutateProteinModelDomain() :
      m_Collector(),
      m_Locator(),
      m_Mutate(),
      m_Scheme( GetStaticClassName< MutateProteinModelDomain>())
    {
    }

    //! @brief constructor from a CollectorInterface and a MutateInterface
    //! @param COLLECTOR Collector that returns Domains from a given ProteinModel
    //! @param MOVE Move that works on a Domain
    //! @param SCHEME Scheme to be used
    MutateProteinModelDomain::MutateProteinModelDomain
    (
      const find::CollectorInterface< util::ShPtrVector< assemble::Domain>, assemble::ProteinModel> &COLLECTOR,
      const coord::MoveInterface &MOVE,
      const std::string &SCHEME
    ) :
      m_Collector( COLLECTOR.Clone()),
      m_Locator(),
      m_Mutate( new math::MutateMoveWrapper< assemble::Domain>( MOVE)),
      m_Scheme( SCHEME)
    {
    }

    //! @brief constructor from a CollectorInterface and a MutateInterface
    //! @param COLLECTOR Collector that returns Domains from a given ProteinModel
    //! @param MUTATE Mutate that works on a Domain
    //! @param SCHEME Scheme to be used
    MutateProteinModelDomain::MutateProteinModelDomain
    (
      const find::CollectorInterface< util::ShPtrVector< assemble::Domain>, assemble::ProteinModel> &COLLECTOR,
      const math::MutateInterface< assemble::Domain> &MUTATE,
      const std::string &SCHEME
    ) :
      m_Collector( COLLECTOR.Clone()),
      m_Locator(),
      m_Mutate( MUTATE.Clone()),
      m_Scheme( SCHEME)
    {
    }

    //! @brief constructor from a LocatorInterface and a MutateInterface
    //! @param LOCATOR Locator that returns a Domain from a given ProteinModel
    //! @param MOVE Move that works on a Domain
    //! @param SCHEME Scheme to be used
    MutateProteinModelDomain::MutateProteinModelDomain
    (
      const find::LocatorInterface< util::ShPtr< assemble::Domain>, assemble::ProteinModel> &LOCATOR,
      const coord::MoveInterface &MOVE,
      const std::string &SCHEME
    ) :
      m_Collector(),
      m_Locator( LOCATOR.Clone()),
      m_Mutate( new math::MutateMoveWrapper< assemble::Domain>( MOVE)),
      m_Scheme( SCHEME)
    {
    }

    //! @brief constructor from a LocatorInterface and a MutateInterface
    //! @param LOCATOR Locator that returns a Domain from a given ProteinModel
    //! @param MUTATE Mutate that works on a Domain
    //! @param SCHEME Scheme to be used
    MutateProteinModelDomain::MutateProteinModelDomain
    (
      const find::LocatorInterface< util::ShPtr< assemble::Domain>, assemble::ProteinModel> &LOCATOR,
      const math::MutateInterface< assemble::Domain> &MUTATE,
      const std::string &SCHEME
    ) :
      m_Collector(),
      m_Locator( LOCATOR.Clone()),
      m_Mutate( MUTATE.Clone()),
      m_Scheme( SCHEME)
    {
    }

    //! @brief Clone function
    //! @return pointer to new MutateProteinModelDomain
    MutateProteinModelDomain *MutateProteinModelDomain::Clone() const
    {
      return new MutateProteinModelDomain( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &MutateProteinModelDomain::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief operator that tkes a ProteinModel and return a mutated ProteinModel
    //! @param PROTEIN_MODEL ProteinModel which will be mutated
    //! @return MutateResult with the mutated ProteinModel
    math::MutateResult< assemble::ProteinModel> MutateProteinModelDomain::operator()
    (
      const assemble::ProteinModel &PROTEIN_MODEL
    ) const
    {
      // static empty result
      static util::ShPtr< assemble::ProteinModel> s_empty_model;

      BCL_MessageVrb( " mutate domain: " + m_Scheme);

      // store ShPtr to domain
      util::ShPtr< assemble::Domain> this_domain;

      // if locator is given
      if( m_Locator.IsDefined())
      {
        this_domain = m_Locator->Locate( PROTEIN_MODEL);
      }
      // if collector is defined
      else if( m_Collector.IsDefined())
      {
        // use the collector to get the domains
        util::ShPtrVector< assemble::Domain> domains( m_Collector->Collect( PROTEIN_MODEL));

        // if there are no domains
        if( domains.IsEmpty())
        {
          BCL_MessageVrb( "No domains found returning failure");
          return math::MutateResult< assemble::ProteinModel>( s_empty_model, *this);
        }
        // pick one of the domains randomly
        this_domain = ( domains( random::GetGlobalRandom().Random< size_t>( domains.GetSize() - 1)));
      }

      // if domain is not defined
      if( !this_domain.IsDefined())
      {
        BCL_MessageVrb( "Found domain is not defined");
        return math::MutateResult< assemble::ProteinModel>( s_empty_model, *this);
      }

      // apply the mutate to this domain and store the result
      math::MutateResult< assemble::Domain> domain_mutate_result( m_Mutate->operator ()( *this_domain));

      // if the mutate was not successful
      if( !domain_mutate_result.GetArgument().IsDefined())
      {
        return math::MutateResult< assemble::ProteinModel>( s_empty_model, *this);
      }
      // otherwise create a new model
      util::ShPtr< assemble::ProteinModel> new_model( PROTEIN_MODEL.Clone());

      // replace the mutated domain within the new_model
      new_model->Replace( domain_mutate_result.GetArgument()->GetData());

      // end
      return math::MutateResult< assemble::ProteinModel>( new_model, *this);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &MutateProteinModelDomain::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Collector, ISTREAM);
      io::Serialize::Read( m_Locator, ISTREAM);
      io::Serialize::Read( m_Mutate, ISTREAM);
      io::Serialize::Read( m_Scheme, ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &MutateProteinModelDomain::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Collector, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Locator, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Mutate, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Scheme, OSTREAM, INDENT) << '\n';

      // end
      return OSTREAM;
    }

  } // namespace fold
} // namespace bcl
