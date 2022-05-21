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

#ifndef BCL_FOLD_MUTATE_PROTEIN_MODEL_DOMAIN_H_
#define BCL_FOLD_MUTATE_PROTEIN_MODEL_DOMAIN_H_

// include the namespace header
#include "bcl_fold.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"
#include "coord/bcl_coord.fwd.hh"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_ss_types.h"
#include "find/bcl_find_collector_interface.h"
#include "find/bcl_find_locator_interface.h"
#include "math/bcl_math_mutate_interface.h"
#include "util/bcl_util_sh_ptr_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MutateProteinModelDomain
    //! @brief class is a wrapper class that allows use of Mutate classes that work on Domains
    //! @details This class allows a modular wrapper for using collectors that return Domains and mutates that work on Domains.
    //! By switching in different versions of collectors and mutates, different copies of this class can be created
    //! that are tailored for different mutate strategies.
    //!
    //! @see @link example_fold_mutate_protein_model_domain.cpp @endlink
    //! @author karakam
    //! @date Mar 12, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MutateProteinModelDomain :
      public math::MutateInterface< assemble::ProteinModel>
    {

    private:

    //////////
    // data //
    //////////

      //! collector that returns a vector of domains from a protein model
      util::ShPtr< find::CollectorInterface< util::ShPtrVector< assemble::Domain>, assemble::ProteinModel> > m_Collector;

      //! locator that returns a domain from a protein model
      util::ShPtr< find::LocatorInterface< util::ShPtr< assemble::Domain>, assemble::ProteinModel> > m_Locator;

      //! move that will be applied
      util::ShPtr< math::MutateInterface< assemble::Domain> > m_Mutate;

      //! scheme
      std::string m_Scheme;

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      MutateProteinModelDomain();

      //! @brief constructor from a CollectorInterface and a MutateInterface
      //! @param COLLECTOR Collector that returns Domains from a given ProteinModel
      //! @param MOVE Move that works on a Domain
      //! @param SCHEME Scheme to be used
      MutateProteinModelDomain
      (
        const find::CollectorInterface< util::ShPtrVector< assemble::Domain>, assemble::ProteinModel> &COLLECTOR,
        const coord::MoveInterface &MOVE,
        const std::string &SCHEME = GetStaticClassName< MutateProteinModelDomain>()
      );

      //! @brief constructor from a CollectorInterface and a MutateInterface
      //! @param COLLECTOR Collector that returns Domains from a given ProteinModel
      //! @param MUTATE Mutate that works on a Domain
      //! @param SCHEME Scheme to be used
      MutateProteinModelDomain
      (
        const find::CollectorInterface< util::ShPtrVector< assemble::Domain>, assemble::ProteinModel> &COLLECTOR,
        const math::MutateInterface< assemble::Domain> &MUTATE,
        const std::string &SCHEME = GetStaticClassName< MutateProteinModelDomain>()
      );

      //! @brief constructor from a LocatorInterface and a MutateInterface
      //! @param LOCATOR Locator that returns a Domain from a given ProteinModel
      //! @param MOVE Move that works on a Domain
      //! @param SCHEME Scheme to be used
      MutateProteinModelDomain
      (
        const find::LocatorInterface< util::ShPtr< assemble::Domain>, assemble::ProteinModel> &LOCATOR,
        const coord::MoveInterface &MOVE,
        const std::string &SCHEME = GetStaticClassName< MutateProteinModelDomain>()
      );

      //! @brief constructor from a LocatorInterface and a MutateInterface
      //! @param LOCATOR Locator that returns a Domain from a given ProteinModel
      //! @param MUTATE Mutate that works on a Domain
      //! @param SCHEME Scheme to be used
      MutateProteinModelDomain
      (
        const find::LocatorInterface< util::ShPtr< assemble::Domain>, assemble::ProteinModel> &LOCATOR,
        const math::MutateInterface< assemble::Domain> &MUTATE,
        const std::string &SCHEME = GetStaticClassName< MutateProteinModelDomain>()
      );

      //! @brief Clone function
      //! @return pointer to new MutateProteinModelDomain
      MutateProteinModelDomain *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief gets the scheme for this mutate
      //! @return the scheme for this mutate
      const std::string &GetScheme() const
      {
        return m_Scheme;
      }

    ///////////////
    // operators //
    ///////////////

      //! @brief operator that takes a ProteinModel and return a mutated ProteinModel
      //! @param PROTEIN_MODEL ProteinModel which will be mutated
      //! @return MutateResult with the mutated ProteinModel
      math::MutateResult< assemble::ProteinModel> operator()( const assemble::ProteinModel &PROTEIN_MODEL) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    }; // class MutateProteinModelDomain

  } // namespace fold
} // namespace bcl

#endif // BCL_FOLD_MUTATE_PROTEIN_MODEL_DOMAIN_H_ 
