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

#ifndef BCL_FOLD_COLLECTOR_LOOP_DOMAIN_H_
#define BCL_FOLD_COLLECTOR_LOOP_DOMAIN_H_

// include the namespace header
#include "bcl_fold.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_atom_types.h"
#include "find/bcl_find_collector_interface.h"
#include "storage/bcl_storage_set.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class CollectorLoopDomain
    //! @brief Collector class for collecting all or only unclosed loop domains from a given protein model
    //! @details This class uses ProteinModelData of the given protein model to access the loop domain locators. Then
    //! it uses these loop domain locators, to locate the loop domains within the given model. The default behaviour,
    //! accessed by using the default constructor or setting the m_CollectUnclosedOnly to false, collects all the
    //! loop domains. If m_CollectUnclosedOnly is set to true, then it will collect only the unclosed loops in the model
    //! where it will use the given m_LoopClosureThreshold and the m_SuperimposeAtomTypes to find out whether each
    //! loop is closed or not.
    //!
    //! @see @link example_fold_collector_loop_domain.cpp @endlink
    //! @author karakam
    //! @date Feb 15, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API CollectorLoopDomain :
      public find::CollectorInterface< util::ShPtrList< LoopDomain>, assemble::DomainInterface>
    {

    private:

    //////////
    // data //
    //////////

      //! boolean to collect only the unclosed ones
      bool m_CollectUnclosedOnly;

      //! loop closure threshold
      double m_LoopClosureThreshold;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      CollectorLoopDomain();

      //! @brief constructor from variables
      //! @param COLLECT_UNCLOSED_ONLY boolean to collect only the unclosed ones
      //! @param LOOP_CLOSURE_THRESHOLD Distance sum threshold for identifying closed loops
      CollectorLoopDomain
      (
        const bool COLLECT_UNCLOSED_ONLY,
        const double LOOP_CLOSURE_THRESHOLD
      );

      //! @brief Clone function
      //! @return pointer to new CollectorLoopDomain
      CollectorLoopDomain *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief get a short name for this class
      //! @return a short name for this class
      const std::string &GetAlias() const;

      //! @brief returns parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief Collect loop domains from given protein model
      //! @param PROTEIN_MODEL ProteinModel of interest
      //! @return loop domains collected from given protein model
      util::ShPtrList< LoopDomain> Collect( const assemble::DomainInterface &PROTEIN_MODEL) const;

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

    }; // class CollectorLoopDomain

  } // namespace fold
} // namespace bcl

#endif // BCL_FOLD_COLLECTOR_LOOP_DOMAIN_H_ 
