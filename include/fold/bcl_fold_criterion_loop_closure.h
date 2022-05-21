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

#ifndef BCL_FOLD_CRITERION_LOOP_CLOSURE_H_
#define BCL_FOLD_CRITERION_LOOP_CLOSURE_H_

// include the namespace header
#include "bcl_fold.h"

// includes from bcl - sorted alphabetically
#include "bcl_fold_locator_loop_domain.h"
#include "command/bcl_command_parameter.h"
#include "command/bcl_command_parameter_check_ranged.h"
#include "command/bcl_command_parameter_interface.h"
#include "opti/bcl_opti_criterion_interface.h"
#include "util/bcl_util_sh_ptr_list.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class CriterionLoopClosure
    //! @brief is a terminate interface derived class for determining if a loop is closed, and terminating if closed
    //! @details It determines loop close by comparing the positions of a set of
    //! atoms in a pseudo residue and the actual residue it corresponds to. It calculates the sum of the squares of the
    //! distances between the pseudo residue atoms and the actual residue atoms. This sum distance is then compared
    //! to a closure threshold. If the sum distance is less than the closure threshold then the terminate is considered
    //! to have its criteria met.
    //!
    //! @tparam t_ArgumentType is the type of the approximation argument
    //! @tparam t_ResultType is the type of the approximation result
    //!
    //! @see @link example_fold_criterion_loop_closure.cpp @endlink
    //! @author alexanns, fischea
    //! @date Feb 12, 2013
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_ArgumentType, typename t_ResultType>
    class CriterionLoopClosure :
      public opti::CriterionInterface< t_ArgumentType, t_ResultType>
    {

    //////////
    // data //
    //////////

    private:

      //! the threshold for considering the loop closed
      double m_RMSDClosureThreshold;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // @brief default constructor
      CriterionLoopClosure() :
        m_RMSDClosureThreshold( 0.0)
      {
      }

      //! @brief construct from loop closure threshold
      //! @param RMSD_CLOSURE_THRESHOLD threshold for considering the loop closed
      CriterionLoopClosure( const double &RMSD_CLOSURE_THRESHOLD) :
        m_RMSDClosureThreshold( RMSD_CLOSURE_THRESHOLD)
      {
      }

      //! @brief Clone function
      //! @return pointer to a new CriterionLoopClosure< t_ArgumentType, t_ResultType>
      CriterionLoopClosure *Clone() const
      {
        return new CriterionLoopClosure< t_ArgumentType, t_ResultType>( *this);
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief returns the class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      //! @brief returns the class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetAlias() const
      {
        static const std::string s_alias( "LoopClosure");
        return s_alias;
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief returns whether the maximum number of iterations has been observed for the given tracker
      //! @param TRACKER tracker to evaluate the criterion for
      //! @return true, if any of the termination criteria are met yet
      bool CriteriaMet( const opti::Tracker< t_ArgumentType, t_ResultType> &TRACKER) const
      {
        // reference to model
        const assemble::ProteinModel &model( TRACKER.GetCurrent()->First());

        // get the loop domain locators from the last protein model in the tracker
        const util::ShPtr< util::ShPtrList< LocatorLoopDomain> > sp_loop_domain_locators
        (
          model.GetProteinModelData()->GetData( assemble::ProteinModelData::e_LoopDomainLocators)
        );

        // make sure the loop domain locators are initialized
        BCL_Assert( sp_loop_domain_locators.IsDefined(), "Loop domain locators are not stored with the protein model");

        // iterate over the loop domain locators
        for
        (
          util::ShPtrList< LocatorLoopDomain>::const_iterator loop_itr( sp_loop_domain_locators->Begin()),
          loop_itr_end( sp_loop_domain_locators->End()); loop_itr != loop_itr_end; ++loop_itr
        )
        {
          if( !LocatorLoopDomain::IsClosed( **loop_itr, model, m_RMSDClosureThreshold))
          {
            return false;
          }
        }

        // if this point is reached, it means all loop ends were succesfully superimposed below the given threshold
        // so the criteria has been met
        return true;
      }

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM)
      {
        // read members
        io::Serialize::Read( m_RMSDClosureThreshold, ISTREAM);

        // end
        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        // write members
        io::Serialize::Write( m_RMSDClosureThreshold, OSTREAM, INDENT);

        // end
        return OSTREAM;
      }

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const
      {
        io::Serializer serializer;
        serializer.SetClassDescription( "Triggers if all loops are closed to within the given rmsd");
        serializer.AddInitializer
        (
          "",
          "RMSD threshold for a loop to be considered closed",
          io::Serialization::GetAgent( &m_RMSDClosureThreshold),
          "0.08"
        );
        return serializer;
      }

    }; // template class CriterionLoopClosure< t_ArgumentType, t_ResultType>

    // instantiate s_Instance
    template< typename t_ArgumentType, typename t_ResultType>
    const util::SiPtr< const util::ObjectInterface> CriterionLoopClosure< t_ArgumentType, t_ResultType>::s_Instance
    (
      util::Enumerated< opti::CriterionInterface< t_ArgumentType, t_ResultType> >::AddInstance( new CriterionLoopClosure< t_ArgumentType, t_ResultType>())
    );

  } // namespace fold
} // namespace bcl

#endif // BCL_FOLD_CRITERION_LOOP_CLOSURE_H_
