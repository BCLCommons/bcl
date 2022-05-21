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

// include example header
#include "example.h"
// include the header of the class which this example is for
#include "signal/bcl_signal_slots.h"

// includes from bcl - sorted alphabetically
#include "signal/bcl_signal_signal.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_signal_slots.cpp
  //!
  //! @author woetzen
  //! @date Nov 6, 2010 
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on 
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleSignalSlots :
    public ExampleInterface
  {
    class Switch
    {

    public:

      signal::Signal0 m_Clicked;

    };

    class Light :
      public signal::Slots
    {
    //////////
    // data //
    //////////

      bool m_State; // true for on, false is off
      std::string m_Name; //! name of light

    public:

      Light( const std::string &NAME) :
        m_State( false),
        m_Name( NAME)
      {}

      void ToggleState()
      {
        m_State = !m_State;
        MessageState();
      }

      void TurnOn()
      {
        m_State = true;
        MessageState();
      }

      void TurnOff()
      {
        m_State = false;
        MessageState();
      }

    private:

      void MessageState() const
      {
        BCL_MessageVrb( m_Name + "light is " + ( m_State ? "on" : "off"));
      }

    }; // class Light

  public:

    ExampleSignalSlots *Clone() const
    {
      return new ExampleSignalSlots( *this);
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

    int Run() const
    {
    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      Light light_kitchen( "kitchen");
      Light light_office( "office");
      Switch switch_kitchen;
      Switch switch_office;
      Switch switch_appartment;

    /////////////////
    // data access //
    /////////////////

      switch_kitchen.m_Clicked.Connect( &light_kitchen, &Light::ToggleState);
      switch_office.m_Clicked.Connect( &light_office, &Light::ToggleState);
      switch_appartment.m_Clicked.Connect( &light_kitchen, &Light::ToggleState);
      switch_appartment.m_Clicked.Connect( &light_office, &Light::ToggleState);

    ///////////////
    // operators //
    ///////////////

      switch_kitchen.m_Clicked.Emit();
      switch_office.m_Clicked.Emit();
      switch_kitchen.m_Clicked.Emit();
      switch_office.m_Clicked.Emit();
      switch_appartment.m_Clicked.Emit();

    ////////////////
    // operations //
    ////////////////

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleSignalSlots

  const ExampleClass::EnumType ExampleSignalSlots::s_Instance
  (
    GetExamples().AddEnum( ExampleSignalSlots())
  );

} // namespace bcl
