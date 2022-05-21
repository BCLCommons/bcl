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

#ifndef BCL_CHEMISTRY_REACTION_ENSEMBLE_H_
#define BCL_CHEMISTRY_REACTION_ENSEMBLE_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_chemistry_reaction_complete.h"
#include "math/bcl_math_range.h"
#include "storage/bcl_storage_list.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ReactionEnsemble
    //! @brief Container class for Reaction objects
    //!
    //! @see @link example_chemistry_reaction_ensemble.cpp @endlink
    //! @author geanesar
    //! @date Jul 30, 2015
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ReactionEnsemble :
      //public storage::List< ReactionComplete>
      public util::ObjectInterface
    {

    private:

    //////////
    // data //
    //////////

      //! list of share pointers to SmallMoleculeConformations
      storage::List< ReactionComplete> m_ReactionEnsemble;

    public:

    //////////////
    // typedefs //
    //////////////

      typedef storage::List< ReactionComplete>::const_iterator const_iterator;
      typedef storage::List< ReactionComplete>::iterator iterator;

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      ReactionEnsemble()
      {}

      //! @brief construct from list of reactions 
      ReactionEnsemble( const storage::List< ReactionComplete> &REACTIONS) :
        m_ReactionEnsemble( REACTIONS)
      {}

      //! @brief construct from an input stream linked to RXN format
      //! @param ISTREAM input stream, reads in RXN format
      //! @param RANGE a range of reactions to load, by default loads all reactions in the stream
      explicit ReactionEnsemble
      (
        std::istream &ISTREAM,
        const math::Range< size_t> &RANGE = math::Range< size_t>( 0, std::numeric_limits< size_t>::max())
      );

      //! @brief Clone function
      //! @return pointer to new ReactionEnsemble
      ReactionEnsemble *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief return a shared pointer vector containing the ensemble
      const storage::List< ReactionComplete> &GetReactions() const
      {
        return m_ReactionEnsemble;
      }

      //! @brief return a shared pointer vector containing the ensemble
      storage::List< ReactionComplete> &GetReactions()
      {
        return m_ReactionEnsemble;
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief add a ReactionComplete to the ensemble
      //! @param RXN_OBJ the object to add to the ensemble
      void PushBack( const ReactionComplete &RXN_OBJ)
      {
        m_ReactionEnsemble.PushBack( RXN_OBJ);
      }

      //! @brief merge this ensemble with another ensemble
      //! @param ENSEMBLE the object to add to the ensemble
      void Append( const ReactionEnsemble &ENSEMBLE)
      {
        m_ReactionEnsemble.Append( ENSEMBLE.m_ReactionEnsemble);
      }

      //! @brief read additional molecules into the ensemble
      //! @param ISTREAM input stream, reads in RXN format
      //! @param RANGE a range of reactions to load, by default loads all reactions in the stream
      void ReadMoreFromRXN
      (
        std::istream &ISTREAM,
        const math::Range< size_t> &RANGE = math::Range< size_t>( 0, std::numeric_limits< size_t>::max())
      );

      //! @brief GetSize returns size of the container
      size_t GetSize() const
      {
        return m_ReactionEnsemble.GetSize();
      }

      //! @brief return true if the ensemble is empty
      bool IsEmpty() const
      {
        return m_ReactionEnsemble.IsEmpty();
      }

      //! @brief const iterator begin
      const_iterator Begin() const
      {
        return m_ReactionEnsemble.Begin();
      }

      //! @brief const iterator end
      const_iterator End() const
      {
        return m_ReactionEnsemble.End();
      }

      //! @brief iterator begin
      iterator Begin()
      {
        return m_ReactionEnsemble.Begin();
      }

      //! @brief iterator end
      iterator End()
      {
        return m_ReactionEnsemble.End();
      }

    ///////////////
    // operators //
    ///////////////

    //////////////////////
    // input and output //
    //////////////////////

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @return output stream which was written to
      std::ostream &WriteRXN( std::ostream &OSTREAM) const;

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      virtual std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return output stream which was written to
      virtual std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    //////////////////////
    // helper functions //
    //////////////////////

    }; // class ReactionEnsemble

  } // namespace chemistry
} // namespace bcl

#endif // BCL_CHEMISTRY_REACTION_ENSEMBLE_H_
