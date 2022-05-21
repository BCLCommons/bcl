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

#ifndef BCL_OPTI_APPROXIMATOR_POWELL_H_
#define BCL_OPTI_APPROXIMATOR_POWELL_H_

// include the namespace header
#include "bcl_opti.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_opti_approximator_golden_section.h"
#include "bcl_opti_criterion_combine.h"
#include "bcl_opti_criterion_convergence_result.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace opti
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ApproximatorPowell
    //! @brief this class implements Powell's method for finding a local minimum
    //! @details an n-dimensional function is optimized using n starting directions, doing line searches in each
    //! direction and applying the combination of each line searche's result and replacing the "worst" direction with
    //! this direction
    //!
    //! @see http://math.fullerton.edu/mathews/n2003/PowellMethodMod.html
    //! @see http://en.wikipedia.org/wiki/Powell%27s_method
    //!
    //! @tparam t_ArgumentType a subtraction has to be defined for that type as well as a multiplication
    //! @tparam t_ResultType a subtraction has to be defined for that type as well as devision by t_ArgumentType and a
    //! multiplication
    //!
    //! @see @link example_opti_approximator_powell.cpp @endlink
    //! @author woetzen, fischea
    //! @date Dec 18, 2012
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_ArgumentType, typename t_ResultType>
    class ApproximatorPowell :
      public ApproximatorModularBase< t_ArgumentType, t_ResultType>
    {

    //////////
    // data //
    //////////

    private:

      //! ShPtr to the objective function
      util::ShPtr< math::FunctionInterfaceSerializable< t_ArgumentType, t_ResultType> > m_Objective;

      //! search directions
      storage::Vector< t_ArgumentType> m_SearchDirections;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      ApproximatorPowell() :
        ApproximatorModularBase< t_ArgumentType, t_ResultType>( e_SmallerAbsIsBetter),
        m_Objective(),
        m_SearchDirections()
      {
      }

      //! @brief construct from members
      //! @param OBJECTIVE the objective function
      //! @param CRITERION the termination criterion
      //! @param ARGUMENT initial argument
      //! @param SEARCH_DIRECTIONS search directions
      ApproximatorPowell
      (
        const math::FunctionInterfaceSerializable< t_ArgumentType, t_ResultType> &OBJECTIVE,
        const CriterionInterface< t_ArgumentType, t_ResultType> &CRITERION,
        const storage::Vector< t_ArgumentType> &SEARCH_DIRECTIONS,
        const t_ArgumentType &ARGUMENT
      ) :
        ApproximatorModularBase< t_ArgumentType, t_ResultType>( e_SmallerAbsIsBetter),
        m_Objective( util::CloneToShPtr( OBJECTIVE)),
        m_SearchDirections( SEARCH_DIRECTIONS)
      {
        // set the termination criterion
        this->SetCriterion( CRITERION);

        // track initial approximation result
        this->Track
        (
          util::ShPtr< storage::Pair< t_ArgumentType, t_ResultType> >
          (
            new storage::Pair< t_ArgumentType, t_ResultType>
            (
              ARGUMENT, m_Objective->operator()( ARGUMENT)
            )
          )
        );
      }

      //! @brief Clone function
      //! @return pointer to a new ApproximatorPowell< t_ArgumentType, t_ResultType>
      ApproximatorPowell *Clone() const
      {
        return new ApproximatorPowell( *this);
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief returns the class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief returns the current approximation
      //! @return current argument result pair
      const util::ShPtr< storage::Pair< t_ArgumentType, t_ResultType> > GetCurrentApproximation() const
      {
        return this->GetTracker().GetCurrent();
      }

    protected:

       //! @brief conducts the next approximation step and stores the approximation
       void Next()
       {
         // get start, previous and current approximation steps
         const util::ShPtr< storage::Pair< t_ArgumentType, t_ResultType> > sp_start( this->GetTracker().GetCurrent());
         util::ShPtr< storage::Pair< t_ArgumentType, t_ResultType> > sp_previous( sp_start);
         util::ShPtr< storage::Pair< t_ArgumentType, t_ResultType> > sp_current;

         // the best direction and the function difference
         typename storage::Vector< t_ArgumentType>::iterator itr_best( m_SearchDirections.Begin());
         t_ResultType best_difference( 0);

         // iterate over all search directions
         for
         (
           typename storage::Vector< t_ArgumentType>::iterator itr
           (
             m_SearchDirections.Begin()), itr_end( m_SearchDirections.End()
           );
           itr != itr_end;
           ++itr, sp_previous = sp_current
         )
         {

           /////////////////////////////////
           // conduct a line search along //
           // the current direction       //
           /////////////////////////////////

           // create termination criteria
           CriterionCombine< t_ArgumentType, t_ResultType> criterion_combine;
           criterion_combine.InsertCriteria( CriterionNumberIterations< t_ArgumentType, t_ResultType>( 15));
           criterion_combine.InsertCriteria( CriterionConvergenceResult< t_ArgumentType, t_ResultType>( 2, 0.0001));

           // define borders for the line search
           const t_ArgumentType right_border( sp_previous->First() - ( *itr));
           const t_ArgumentType left_border( sp_previous->First() + ( *itr));

           // create approximator for the line search
           ApproximatorGoldenSection< t_ArgumentType, t_ResultType> line_search_approximator
           (
             *m_Objective,
             criterion_combine,
             left_border,
             right_border
           );

           // approximate
           line_search_approximator.Approximate();

           // update current
           sp_current = line_search_approximator.GetTracker().GetBest();

           // line search improved previous result
           if( sp_current->Second() < sp_previous->Second())
           {
             // insert the new result into the tracker
             this->Track( sp_current);
           }
           else // line search did not help
           {
             BCL_Message( util::Message::e_Verbose, "no improvement in line search");
             sp_current = sp_previous;

             // shorten the direction
             ( *itr) *= 0.5;
           }

          // difference
          const t_ResultType difference( sp_previous->Second() - sp_current->Second());
          if( difference > best_difference)
          {
            itr_best = itr;
            best_difference = difference;
          }

         } // line search in all directions

         // extended direction
         util::ShPtr< storage::Pair< t_ArgumentType, t_ResultType> > sp_extended
         (
           new storage::Pair< t_ArgumentType, t_ResultType>
           (
             2.0 * sp_current->First() - sp_start->First(), t_ResultType()
           )
         );
         sp_extended->Second() = m_Objective->operator()( sp_extended->First());

         // if there is no further decrease in the average direction expected, no direction updates
         if( sp_extended->Second() >= sp_start->Second())
         {
           this->Track( sp_current);
         }
         else if
         (
           2 * ( sp_start->Second() - 2 * sp_current->Second() + sp_extended->Second()) *
           math::Sqr( sp_start->Second() - sp_current->Second() - best_difference) >=
           ( best_difference * math::Sqr( sp_start->Second() - sp_extended->Second()))
         )
         {
           this->Track( sp_current);
         }
         else
         {
           // replace the best search direction with overall direction
           BCL_Message( util::Message::e_Verbose, "update the best direction with the average direction");
           *itr_best = sp_current->First() - sp_start->First();

           // place the average direction first
           if( itr_best != m_SearchDirections.Begin())
           {
             std::swap( *itr_best, m_SearchDirections.FirstElement());
           }

           // track pair of argument and function value
           this->Track( sp_current);
         }

         BCL_Message( util::Message::e_Standard, "CURRENT: " + util::Format()( *sp_current));
       }

       //! @brief evaluates whether the approximation can continue
       //! @return true, if the approximation can continue - otherwise false
       bool CanContinue() const
       {
         return true;
       }

    //////////////////////
    // input and output //
    //////////////////////

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM)
      {
        // call Read of base class
        ApproximatorModularBase< t_ArgumentType, t_ResultType>::Read( ISTREAM);

        // read members
        io::Serialize::Read( m_Objective, ISTREAM);
        io::Serialize::Read( m_SearchDirections, ISTREAM);

        // end
        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        // call Write of base class
        ApproximatorModularBase< t_ArgumentType, t_ResultType>::Write( OSTREAM, INDENT) << '\n';

        // write members
        io::Serialize::Write( m_Objective, OSTREAM, INDENT) << '\n';
        io::Serialize::Write( m_SearchDirections, OSTREAM, INDENT);

        // end
        return OSTREAM;
      }

    }; // template class ApproximatorPowell< t_ArgumentType, t_ResultType>

  } // namespace opti
} // namespace bcl

#endif // BCL_OPTI_APPROXIMATOR_POWELL_H_
