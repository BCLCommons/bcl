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

#ifndef BCL_NMR_SIGNAL_H_
#define BCL_NMR_SIGNAL_H_

// include the namespace header
#include "bcl_nmr.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_nmr_signal_1d.h"
#include "util/bcl_util_sh_ptr_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace nmr
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Signal
    //! @brief consisting of 1D and ND signals are constructed with using bcl_signal.
    //!
    //! @see @link example_nmr_signal.cpp @endlink
    //! @author mueller, butkiem1
    //! @date 2006/05/16
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API Signal :
      public util::ObjectInterface
    {
    private:

    //////////
    // data //
    //////////

      util::ShPtrVector< Signal1D> m_Signals1D; //!< stores the signal for every dimension of the spectra

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief standard constructor
      Signal()
      {
      }

      //! @brief copy constructor from other signal objects properties
      Signal( const Signal &SIGNAL) :
        m_Signals1D( SIGNAL.m_Signals1D)
      {
      }

      //! @brief construct signal from its properties
      Signal( const util::ShPtrVector< Signal1D> &SIGNALS1D) :
        m_Signals1D( SIGNALS1D)
      {
      }

      //! @brief virtual copy constructor
      Signal *Clone() const
      {
        return new Signal( *this);
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      //! @brief get the stored 1D signals
      //! @return ShPtrVector on the stored 1D signals
      const util::ShPtrVector< Signal1D> &GetSignals1D() const
      {
        return m_Signals1D;
      }

      //! @brief set 1D signals
      void SetSignals1D( const util::ShPtrVector< Signal1D> &SIGNALS1D)
      {
        m_Signals1D = SIGNALS1D;
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief check if given atom is involved in that signal
      //! @param ATOM the atom of interest
      //! @return true, if any Signal1D has the given atom
      bool ContainsAtom( const chemistry::AtomConformationalInterface &ATOM) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! read Spectrum from std::istream
      virtual std::istream &Read( std::istream &ISTREAM);

      //! write Spectrum to std::ostream
      virtual std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    };

  } // namespace nmr
} // namespace bcl

#endif //BCL_NMR_SIGNAL_H_

