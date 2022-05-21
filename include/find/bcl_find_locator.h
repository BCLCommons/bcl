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

#ifndef BCL_FIND_LOCATOR_H_
#define BCL_FIND_LOCATOR_H_

// include the namespace header
#include "bcl_find.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_find_collector_interface.h"
#include "bcl_find_locator_interface.h"
#include "bcl_find_pick_interface.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace find
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Locator
    //! @brief template class to combine a Collector and a Picker
    //! @details This generic Locator class combines a Collector and a Picker
    //!
    //! @see @link example_find_locator.cpp @endlink
    //! @author karakam
    //! @date Feb 21, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_ReturnType, typename t_ArgumentType, typename t_IntermediateResultType>
    class Locator :
      public LocatorInterface< t_ReturnType, t_ArgumentType>
    {

    private:

    //////////
    // data //
    //////////

      //! method for collecting arguments from t_ArgumentType
      util::ShPtr< CollectorInterface< t_IntermediateResultType, t_ArgumentType> > m_Collector;

      //! method for picking argument from the results of m_Collector
      util::ShPtr< PickInterface< t_ReturnType, t_IntermediateResultType> > m_Picker;

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
      Locator() :
        m_Collector(),
        m_Picker()
      {
      }

      //! @brief constructor from a collector and picker
      //! @param COLLECTOR Collector to be used
      //! @param PICKER Picker to be used
      Locator
      (
        const CollectorInterface< t_IntermediateResultType, t_ArgumentType> &COLLECTOR,
        const PickInterface< t_ReturnType, t_IntermediateResultType> &PICKER
      ) :
        m_Collector( COLLECTOR.Clone()),
        m_Picker( PICKER.Clone())
      {
      }

      //! @brief constructor from ShPtrs to a collector and picker
      //! @param SP_COLLECTOR Collector to be used
      //! @param SP_PICKER Picker to be used
      Locator
      (
        const util::ShPtr< CollectorInterface< t_IntermediateResultType, t_ArgumentType> > &SP_COLLECTOR,
        const util::ShPtr< PickInterface< t_ReturnType, t_IntermediateResultType> > &SP_PICKER
      ) :
        m_Collector( SP_COLLECTOR),
        m_Picker( SP_PICKER)
      {
      }

      //! @brief Clone function
      //! @return pointer to new Locator
      virtual Locator *Clone() const
      {
        return new Locator( *this);
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

    ////////////////
    // operations //
    ////////////////

      //! locate the t_ReturnType in t_ArgumentType
      //! @param ARGUMENT entity that contains a t_ReturnType
      //! @return returns the located t_ReturnType
      t_ReturnType Locate( const t_ArgumentType &ARGUMENT) const
      {
        // collect subset and pick one from it
        return m_Picker->Pick( m_Collector->Collect( ARGUMENT));
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
        io::Serialize::Read( m_Collector, ISTREAM);
        io::Serialize::Read( m_Picker, ISTREAM);

        // end
        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        // write members
        io::Serialize::Write( m_Collector, OSTREAM, INDENT) << '\n';
        io::Serialize::Write( m_Picker, OSTREAM, INDENT);

        // end
        return OSTREAM;
      }

    }; // template class Locator

    // instantiate s_Instance
    template
    <
      typename t_ReturnType, typename t_ArgumentType, typename t_IntermediateResultType
    >
    const util::SiPtr< const util::ObjectInterface>
    Locator< t_ReturnType, t_ArgumentType, t_IntermediateResultType>::s_Instance
    (
      GetObjectInstances().AddInstance
      (
        new Locator< t_ReturnType, t_ArgumentType, t_IntermediateResultType>()
      )
    );

  } // namespace find
} // namespace bcl

#endif // BCL_FIND_LOCATOR_H_ 
