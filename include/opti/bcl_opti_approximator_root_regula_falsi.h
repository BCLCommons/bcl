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

#ifndef BCL_OPTI_APPROXIMATOR_ROOT_REGULA_FALSI_H_
#define BCL_OPTI_APPROXIMATOR_ROOT_REGULA_FALSI_H_

// include the namespace header
#include "bcl_opti.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_opti_approximator_modular_base.h"
#include "math/bcl_math_function_interface_serializable.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace opti
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ApproximatorRootRegulaFalsi
    //! @brief Implementation of the RegulaFalsi method to find a root of a continuous function without a first
    //! derivative
    //! @details an objective function and an initial guess for a left and right side argument to the result are
    //! required. this methods unites the advantages of the bisection method and the secant method.
    //!
    //! @see http://en.wikipedia.org/wiki/False_position_method
    //!
    //! @tparam t_ArgumentType a subtraction has to be defined for that type as well as a multiplication
    //! @tparam t_ResultType a subtraction has to be defined for that type as well as devision by t_ArgumentType and a
    //! multiplication
    //!
    //! @see @link example_opti_approximator_root_regula_falsi.cpp @endlink
    //! @author woetzen, fischea
    //! @date Dec 13, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_ArgumentType, typename t_ResultType>
     class ApproximatorRootRegulaFalsi :
       public ApproximatorModularBase< t_ArgumentType, t_ResultType>
     {

    //////////
    // data //
    //////////

    private:

      //! ShPtr to the objective function
      util::ShPtr< math::FunctionInterfaceSerializable< t_ArgumentType, t_ResultType> > m_ObjectiveFunction;

      //! ShPtr to the previous approximation step
      util::ShPtr< storage::Pair< t_ArgumentType, t_ResultType> > m_Previous;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      ApproximatorRootRegulaFalsi() :
        ApproximatorModularBase< t_ArgumentType, t_ResultType>( e_SmallerAbsIsBetter),
        m_ObjectiveFunction(),
        m_Previous()
      {
      }

      //! @brief construct from function, termination criterion and borders of the approximation interval
      //! @param OBJECTIVE_FUNCTION the objective function
      //! @param CRITERION the termination criterion
      //! @param BORDER_LEFT left border of the interval
      //! @param BORDER_RIGHT right border of the interval
      ApproximatorRootRegulaFalsi
      (
        const math::FunctionInterfaceSerializable< t_ArgumentType, t_ResultType> &OBJECTIVE_FUNCTION,
        const CriterionInterface< t_ArgumentType, t_ResultType> &CRITERION,
        const t_ArgumentType &BORDER_LEFT,
        const t_ArgumentType &BORDER_RIGHT
      ) :
        ApproximatorModularBase< t_ArgumentType, t_ResultType>( e_SmallerAbsIsBetter),
        m_ObjectiveFunction( util::CloneToShPtr( OBJECTIVE_FUNCTION)),
        m_Previous
        (
            new storage::Pair< t_ArgumentType, t_ResultType>
            (
              BORDER_LEFT, m_ObjectiveFunction->operator()( BORDER_LEFT)
            )
        )
      {
        // set the termination criterion
        this->SetCriterion( CRITERION);

        // track the initial argument
        this->Track
        (
          util::ShPtr< storage::Pair< t_ArgumentType, t_ResultType> >
          (
            new storage::Pair< t_ArgumentType, t_ResultType>
            (
              BORDER_RIGHT, m_ObjectiveFunction->operator()( BORDER_RIGHT)
            )
          )
        );
      }

      //! @brief Clone function
      //! @return pointer to a new ApproximatorRootRegulaFalsi< t_ArgumentType, t_ResultType>
      ApproximatorRootRegulaFalsi *Clone() const
      {
        return new ApproximatorRootRegulaFalsi( *this);
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

    protected:

      //! @brief conducts the next approximation step and stores the approximation
      void Next()
      {
        // current model
        util::ShPtr< storage::Pair< t_ArgumentType, t_ResultType> > sp_current( this->GetTracker().GetCurrent());

        // bisect given interval
        const t_ArgumentType next_arg
        (
          ( sp_current->First() * m_Previous->Second() - m_Previous->First() * sp_current->Second()) /
          ( m_Previous->Second() - sp_current->Second())
        );

        util::ShPtr< storage::Pair< t_ArgumentType, t_ResultType> > sp_next
        (
          new storage::Pair< t_ArgumentType, t_ResultType>
          (
            next_arg, m_ObjectiveFunction->operator()( next_arg)
          )
        );

        // current and next do have different sign - decrease interval
        if( ( sp_current->Second() < t_ResultType( 0)) != ( sp_next->Second() < t_ResultType( 0)))
        {
          m_Previous = sp_current;
        }
        else
        {
          // correction factor
          const double result_correction( 0.5);

          // correct the previous result
          const t_ResultType corrected_previous_result( result_correction * m_Previous->Second());

          // corrected bisect of given interval
          const t_ArgumentType previous_corrected_arg
          (
            ( sp_current->First() * corrected_previous_result - m_Previous->First() * sp_current->Second()) /
            ( corrected_previous_result - sp_current->Second())
          );

          // pair of argument and function value
          m_Previous = util::ShPtr< storage::Pair< t_ArgumentType, t_ResultType> >
          (
            new storage::Pair< t_ArgumentType, t_ResultType>
            (
              previous_corrected_arg, m_ObjectiveFunction->operator()( previous_corrected_arg)
            )
          );
        }

        // store the new approximation step
        this->Track( sp_next);
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
        io::Serialize::Read( m_ObjectiveFunction, ISTREAM);

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
        io::Serialize::Write( m_ObjectiveFunction, OSTREAM, INDENT) << '\n';

        // end
        return OSTREAM;
      }

    }; // template class ApproximatorRootRegulaFalsi< t_ArgumentType, t_ResultType>

  } // namespace opti
} // namespace bcl

#endif // BCL_OPTI_APPROXIMATOR_ROOT_REGULA_FALSI_H_
