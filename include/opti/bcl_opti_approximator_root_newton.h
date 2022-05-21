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

#ifndef BCL_OPTI_APPROXIMATOR_ROOT_NEWTON_H_
#define BCL_OPTI_APPROXIMATOR_ROOT_NEWTON_H_

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
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ApproximatorRootNewton
    //! @brief Implementation of Newton's Method to find a root for a continuous function with the first derivative
    //!
    //! @tparam t_ArgumentResultType is the type of the approximation argument and result
    //!
    //! @see http://en.wikipedia.org/wiki/Newton%27s_method
    //!
    //! @see @link example_opti_approximator_root_newton.cpp @endlink
    //! @author woetzen, fischea
    //! @date Dec 13, 2012
    //!
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_ArgumentResultType>
    class ApproximatorRootNewton :
      public ApproximatorModularBase< t_ArgumentResultType, t_ArgumentResultType>
    {

    //////////
    // data //
    //////////

    private:

      //! ShPtr to the objective function
      util::ShPtr< math::FunctionInterfaceSerializable< t_ArgumentResultType, t_ArgumentResultType> > m_Objective;

      //! ShPtr to the first derivative of the objective function
      util::ShPtr< math::FunctionInterfaceSerializable< t_ArgumentResultType, t_ArgumentResultType> > m_Derivative;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      ApproximatorRootNewton() :
        ApproximatorModularBase< t_ArgumentResultType, t_ArgumentResultType>( e_SmallerAbsIsBetter),
        m_Objective(),
        m_Derivative()
      {
      }

      //! @brief construct from function, derivative, termination criterion and initial argument
      //! @param OBJECTIVE the objective function
      //! @param DERIVATIVE the first derivative of the objective function
      //! @param CRITERION the termination criterion
      //! @param GUESS initial guess for the root
      ApproximatorRootNewton
      (
        const math::FunctionInterfaceSerializable< t_ArgumentResultType, t_ArgumentResultType> &OBJECTIVE,
        const math::FunctionInterfaceSerializable< t_ArgumentResultType, t_ArgumentResultType> &DERIVATIVE,
        const CriterionInterface< t_ArgumentResultType, t_ArgumentResultType> &CRITERION,
        const t_ArgumentResultType &GUESS
      ) :
        ApproximatorModularBase< t_ArgumentResultType, t_ArgumentResultType>( e_SmallerAbsIsBetter),
        m_Objective( util::CloneToShPtr( OBJECTIVE)),
        m_Derivative( util::CloneToShPtr( DERIVATIVE))
      {
        // set the termination criterion
        this->SetCriterion( CRITERION);

        // track the initial guess
        this->Track
        (
          util::CloneToShPtr
          (
            storage::Pair< t_ArgumentResultType, t_ArgumentResultType>( GUESS, m_Objective->operator()( GUESS))
          )
        );
      }

      //! @brief Clone function
      //! @return pointer to a new ApproximatorRootNewton< t_ArgumentResultType>
      ApproximatorRootNewton *Clone() const
      {
        return new ApproximatorRootNewton( *this);
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
        // current argument and result
        const util::ShPtr< storage::Pair< t_ArgumentResultType, t_ArgumentResultType> > sp_current
        (
          this->GetTracker().GetCurrent()
        );

        // determine next argument Xn+1 = Xn - f(Xn) / f'(Xn)
        util::ShPtr< storage::Pair< t_ArgumentResultType, t_ArgumentResultType> > sp_next_argument_result_value
        (
          new storage::Pair< t_ArgumentResultType, t_ArgumentResultType>
          (
            sp_current->First() - ( sp_current->Second() / m_Derivative->operator()( sp_current->First())),
            t_ArgumentResultType()
          )
        );

        // calculate according next function value
        sp_next_argument_result_value->Second() = m_Objective->operator()( sp_next_argument_result_value->First());

        // store the new approximation step
        this->Track( sp_next_argument_result_value);
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
        ApproximatorModularBase< t_ArgumentResultType, t_ArgumentResultType>::Read( ISTREAM);

        // read members
        io::Serialize::Read( m_Objective, ISTREAM);
        io::Serialize::Read( m_Derivative, ISTREAM);

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
        ApproximatorModularBase< t_ArgumentResultType, t_ArgumentResultType>::Write( OSTREAM, INDENT) << '\n';

        // write members
        io::Serialize::Write( m_Objective, OSTREAM, INDENT) << '\n';
        io::Serialize::Write( m_Derivative, OSTREAM, INDENT);

        // end
        return OSTREAM;
      }

    }; // template class ApproximatorRootNewton< t_ArgumentResultType>

  } // namespace opti
} // namespace bcl

#endif // BCL_OPTI_APPROXIMATOR_ROOT_NEWTON_H_
