// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME TUnfoldV17Dict

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "./TUnfoldIterativeEM.h"
#include "./TUnfoldDensity.h"
#include "./TUnfoldBinningXML.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_TUnfoldV17(void *p = 0);
   static void *newArray_TUnfoldV17(Long_t size, void *p);
   static void delete_TUnfoldV17(void *p);
   static void deleteArray_TUnfoldV17(void *p);
   static void destruct_TUnfoldV17(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TUnfoldV17*)
   {
      ::TUnfoldV17 *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TUnfoldV17 >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TUnfoldV17", ::TUnfoldV17::Class_Version(), "TUnfold.h", 108,
                  typeid(::TUnfoldV17), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TUnfoldV17::Dictionary, isa_proxy, 4,
                  sizeof(::TUnfoldV17) );
      instance.SetNew(&new_TUnfoldV17);
      instance.SetNewArray(&newArray_TUnfoldV17);
      instance.SetDelete(&delete_TUnfoldV17);
      instance.SetDeleteArray(&deleteArray_TUnfoldV17);
      instance.SetDestructor(&destruct_TUnfoldV17);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TUnfoldV17*)
   {
      return GenerateInitInstanceLocal((::TUnfoldV17*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::TUnfoldV17*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_TUnfoldIterativeEMV17(void *p = 0);
   static void *newArray_TUnfoldIterativeEMV17(Long_t size, void *p);
   static void delete_TUnfoldIterativeEMV17(void *p);
   static void deleteArray_TUnfoldIterativeEMV17(void *p);
   static void destruct_TUnfoldIterativeEMV17(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TUnfoldIterativeEMV17*)
   {
      ::TUnfoldIterativeEMV17 *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TUnfoldIterativeEMV17 >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TUnfoldIterativeEMV17", ::TUnfoldIterativeEMV17::Class_Version(), "", 52,
                  typeid(::TUnfoldIterativeEMV17), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TUnfoldIterativeEMV17::Dictionary, isa_proxy, 4,
                  sizeof(::TUnfoldIterativeEMV17) );
      instance.SetNew(&new_TUnfoldIterativeEMV17);
      instance.SetNewArray(&newArray_TUnfoldIterativeEMV17);
      instance.SetDelete(&delete_TUnfoldIterativeEMV17);
      instance.SetDeleteArray(&deleteArray_TUnfoldIterativeEMV17);
      instance.SetDestructor(&destruct_TUnfoldIterativeEMV17);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TUnfoldIterativeEMV17*)
   {
      return GenerateInitInstanceLocal((::TUnfoldIterativeEMV17*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::TUnfoldIterativeEMV17*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_TUnfoldSysV17(void *p = 0);
   static void *newArray_TUnfoldSysV17(Long_t size, void *p);
   static void delete_TUnfoldSysV17(void *p);
   static void deleteArray_TUnfoldSysV17(void *p);
   static void destruct_TUnfoldSysV17(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TUnfoldSysV17*)
   {
      ::TUnfoldSysV17 *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TUnfoldSysV17 >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TUnfoldSysV17", ::TUnfoldSysV17::Class_Version(), "TUnfoldSys.h", 60,
                  typeid(::TUnfoldSysV17), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TUnfoldSysV17::Dictionary, isa_proxy, 4,
                  sizeof(::TUnfoldSysV17) );
      instance.SetNew(&new_TUnfoldSysV17);
      instance.SetNewArray(&newArray_TUnfoldSysV17);
      instance.SetDelete(&delete_TUnfoldSysV17);
      instance.SetDeleteArray(&deleteArray_TUnfoldSysV17);
      instance.SetDestructor(&destruct_TUnfoldSysV17);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TUnfoldSysV17*)
   {
      return GenerateInitInstanceLocal((::TUnfoldSysV17*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::TUnfoldSysV17*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_TUnfoldBinningV17(void *p = 0);
   static void *newArray_TUnfoldBinningV17(Long_t size, void *p);
   static void delete_TUnfoldBinningV17(void *p);
   static void deleteArray_TUnfoldBinningV17(void *p);
   static void destruct_TUnfoldBinningV17(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TUnfoldBinningV17*)
   {
      ::TUnfoldBinningV17 *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TUnfoldBinningV17 >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TUnfoldBinningV17", ::TUnfoldBinningV17::Class_Version(), "TUnfoldBinning.h", 59,
                  typeid(::TUnfoldBinningV17), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TUnfoldBinningV17::Dictionary, isa_proxy, 4,
                  sizeof(::TUnfoldBinningV17) );
      instance.SetNew(&new_TUnfoldBinningV17);
      instance.SetNewArray(&newArray_TUnfoldBinningV17);
      instance.SetDelete(&delete_TUnfoldBinningV17);
      instance.SetDeleteArray(&deleteArray_TUnfoldBinningV17);
      instance.SetDestructor(&destruct_TUnfoldBinningV17);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TUnfoldBinningV17*)
   {
      return GenerateInitInstanceLocal((::TUnfoldBinningV17*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::TUnfoldBinningV17*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_TUnfoldDensityV17(void *p = 0);
   static void *newArray_TUnfoldDensityV17(Long_t size, void *p);
   static void delete_TUnfoldDensityV17(void *p);
   static void deleteArray_TUnfoldDensityV17(void *p);
   static void destruct_TUnfoldDensityV17(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TUnfoldDensityV17*)
   {
      ::TUnfoldDensityV17 *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TUnfoldDensityV17 >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TUnfoldDensityV17", ::TUnfoldDensityV17::Class_Version(), "", 150,
                  typeid(::TUnfoldDensityV17), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TUnfoldDensityV17::Dictionary, isa_proxy, 4,
                  sizeof(::TUnfoldDensityV17) );
      instance.SetNew(&new_TUnfoldDensityV17);
      instance.SetNewArray(&newArray_TUnfoldDensityV17);
      instance.SetDelete(&delete_TUnfoldDensityV17);
      instance.SetDeleteArray(&deleteArray_TUnfoldDensityV17);
      instance.SetDestructor(&destruct_TUnfoldDensityV17);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TUnfoldDensityV17*)
   {
      return GenerateInitInstanceLocal((::TUnfoldDensityV17*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::TUnfoldDensityV17*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_TUnfoldBinningXMLV17(void *p = 0);
   static void *newArray_TUnfoldBinningXMLV17(Long_t size, void *p);
   static void delete_TUnfoldBinningXMLV17(void *p);
   static void deleteArray_TUnfoldBinningXMLV17(void *p);
   static void destruct_TUnfoldBinningXMLV17(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TUnfoldBinningXMLV17*)
   {
      ::TUnfoldBinningXMLV17 *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TUnfoldBinningXMLV17 >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TUnfoldBinningXMLV17", ::TUnfoldBinningXMLV17::Class_Version(), "", 366,
                  typeid(::TUnfoldBinningXMLV17), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TUnfoldBinningXMLV17::Dictionary, isa_proxy, 4,
                  sizeof(::TUnfoldBinningXMLV17) );
      instance.SetNew(&new_TUnfoldBinningXMLV17);
      instance.SetNewArray(&newArray_TUnfoldBinningXMLV17);
      instance.SetDelete(&delete_TUnfoldBinningXMLV17);
      instance.SetDeleteArray(&deleteArray_TUnfoldBinningXMLV17);
      instance.SetDestructor(&destruct_TUnfoldBinningXMLV17);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TUnfoldBinningXMLV17*)
   {
      return GenerateInitInstanceLocal((::TUnfoldBinningXMLV17*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::TUnfoldBinningXMLV17*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr TUnfoldV17::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TUnfoldV17::Class_Name()
{
   return "TUnfoldV17";
}

//______________________________________________________________________________
const char *TUnfoldV17::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TUnfoldV17*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TUnfoldV17::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TUnfoldV17*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TUnfoldV17::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TUnfoldV17*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TUnfoldV17::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TUnfoldV17*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr TUnfoldIterativeEMV17::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TUnfoldIterativeEMV17::Class_Name()
{
   return "TUnfoldIterativeEMV17";
}

//______________________________________________________________________________
const char *TUnfoldIterativeEMV17::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TUnfoldIterativeEMV17*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TUnfoldIterativeEMV17::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TUnfoldIterativeEMV17*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TUnfoldIterativeEMV17::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TUnfoldIterativeEMV17*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TUnfoldIterativeEMV17::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TUnfoldIterativeEMV17*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr TUnfoldSysV17::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TUnfoldSysV17::Class_Name()
{
   return "TUnfoldSysV17";
}

//______________________________________________________________________________
const char *TUnfoldSysV17::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TUnfoldSysV17*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TUnfoldSysV17::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TUnfoldSysV17*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TUnfoldSysV17::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TUnfoldSysV17*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TUnfoldSysV17::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TUnfoldSysV17*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr TUnfoldBinningV17::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TUnfoldBinningV17::Class_Name()
{
   return "TUnfoldBinningV17";
}

//______________________________________________________________________________
const char *TUnfoldBinningV17::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TUnfoldBinningV17*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TUnfoldBinningV17::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TUnfoldBinningV17*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TUnfoldBinningV17::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TUnfoldBinningV17*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TUnfoldBinningV17::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TUnfoldBinningV17*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr TUnfoldDensityV17::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TUnfoldDensityV17::Class_Name()
{
   return "TUnfoldDensityV17";
}

//______________________________________________________________________________
const char *TUnfoldDensityV17::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TUnfoldDensityV17*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TUnfoldDensityV17::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TUnfoldDensityV17*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TUnfoldDensityV17::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TUnfoldDensityV17*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TUnfoldDensityV17::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TUnfoldDensityV17*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr TUnfoldBinningXMLV17::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TUnfoldBinningXMLV17::Class_Name()
{
   return "TUnfoldBinningXMLV17";
}

//______________________________________________________________________________
const char *TUnfoldBinningXMLV17::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TUnfoldBinningXMLV17*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TUnfoldBinningXMLV17::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TUnfoldBinningXMLV17*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TUnfoldBinningXMLV17::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TUnfoldBinningXMLV17*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TUnfoldBinningXMLV17::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TUnfoldBinningXMLV17*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void TUnfoldV17::Streamer(TBuffer &R__b)
{
   // Stream an object of class TUnfoldV17.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TUnfoldV17::Class(),this);
   } else {
      R__b.WriteClassBuffer(TUnfoldV17::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TUnfoldV17(void *p) {
      return  p ? new(p) ::TUnfoldV17 : new ::TUnfoldV17;
   }
   static void *newArray_TUnfoldV17(Long_t nElements, void *p) {
      return p ? new(p) ::TUnfoldV17[nElements] : new ::TUnfoldV17[nElements];
   }
   // Wrapper around operator delete
   static void delete_TUnfoldV17(void *p) {
      delete ((::TUnfoldV17*)p);
   }
   static void deleteArray_TUnfoldV17(void *p) {
      delete [] ((::TUnfoldV17*)p);
   }
   static void destruct_TUnfoldV17(void *p) {
      typedef ::TUnfoldV17 current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::TUnfoldV17

//______________________________________________________________________________
void TUnfoldIterativeEMV17::Streamer(TBuffer &R__b)
{
   // Stream an object of class TUnfoldIterativeEMV17.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TUnfoldIterativeEMV17::Class(),this);
   } else {
      R__b.WriteClassBuffer(TUnfoldIterativeEMV17::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TUnfoldIterativeEMV17(void *p) {
      return  p ? new(p) ::TUnfoldIterativeEMV17 : new ::TUnfoldIterativeEMV17;
   }
   static void *newArray_TUnfoldIterativeEMV17(Long_t nElements, void *p) {
      return p ? new(p) ::TUnfoldIterativeEMV17[nElements] : new ::TUnfoldIterativeEMV17[nElements];
   }
   // Wrapper around operator delete
   static void delete_TUnfoldIterativeEMV17(void *p) {
      delete ((::TUnfoldIterativeEMV17*)p);
   }
   static void deleteArray_TUnfoldIterativeEMV17(void *p) {
      delete [] ((::TUnfoldIterativeEMV17*)p);
   }
   static void destruct_TUnfoldIterativeEMV17(void *p) {
      typedef ::TUnfoldIterativeEMV17 current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::TUnfoldIterativeEMV17

//______________________________________________________________________________
void TUnfoldSysV17::Streamer(TBuffer &R__b)
{
   // Stream an object of class TUnfoldSysV17.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TUnfoldSysV17::Class(),this);
   } else {
      R__b.WriteClassBuffer(TUnfoldSysV17::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TUnfoldSysV17(void *p) {
      return  p ? new(p) ::TUnfoldSysV17 : new ::TUnfoldSysV17;
   }
   static void *newArray_TUnfoldSysV17(Long_t nElements, void *p) {
      return p ? new(p) ::TUnfoldSysV17[nElements] : new ::TUnfoldSysV17[nElements];
   }
   // Wrapper around operator delete
   static void delete_TUnfoldSysV17(void *p) {
      delete ((::TUnfoldSysV17*)p);
   }
   static void deleteArray_TUnfoldSysV17(void *p) {
      delete [] ((::TUnfoldSysV17*)p);
   }
   static void destruct_TUnfoldSysV17(void *p) {
      typedef ::TUnfoldSysV17 current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::TUnfoldSysV17

//______________________________________________________________________________
void TUnfoldBinningV17::Streamer(TBuffer &R__b)
{
   // Stream an object of class TUnfoldBinningV17.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TUnfoldBinningV17::Class(),this);
   } else {
      R__b.WriteClassBuffer(TUnfoldBinningV17::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TUnfoldBinningV17(void *p) {
      return  p ? new(p) ::TUnfoldBinningV17 : new ::TUnfoldBinningV17;
   }
   static void *newArray_TUnfoldBinningV17(Long_t nElements, void *p) {
      return p ? new(p) ::TUnfoldBinningV17[nElements] : new ::TUnfoldBinningV17[nElements];
   }
   // Wrapper around operator delete
   static void delete_TUnfoldBinningV17(void *p) {
      delete ((::TUnfoldBinningV17*)p);
   }
   static void deleteArray_TUnfoldBinningV17(void *p) {
      delete [] ((::TUnfoldBinningV17*)p);
   }
   static void destruct_TUnfoldBinningV17(void *p) {
      typedef ::TUnfoldBinningV17 current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::TUnfoldBinningV17

//______________________________________________________________________________
void TUnfoldDensityV17::Streamer(TBuffer &R__b)
{
   // Stream an object of class TUnfoldDensityV17.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TUnfoldDensityV17::Class(),this);
   } else {
      R__b.WriteClassBuffer(TUnfoldDensityV17::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TUnfoldDensityV17(void *p) {
      return  p ? new(p) ::TUnfoldDensityV17 : new ::TUnfoldDensityV17;
   }
   static void *newArray_TUnfoldDensityV17(Long_t nElements, void *p) {
      return p ? new(p) ::TUnfoldDensityV17[nElements] : new ::TUnfoldDensityV17[nElements];
   }
   // Wrapper around operator delete
   static void delete_TUnfoldDensityV17(void *p) {
      delete ((::TUnfoldDensityV17*)p);
   }
   static void deleteArray_TUnfoldDensityV17(void *p) {
      delete [] ((::TUnfoldDensityV17*)p);
   }
   static void destruct_TUnfoldDensityV17(void *p) {
      typedef ::TUnfoldDensityV17 current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::TUnfoldDensityV17

//______________________________________________________________________________
void TUnfoldBinningXMLV17::Streamer(TBuffer &R__b)
{
   // Stream an object of class TUnfoldBinningXMLV17.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TUnfoldBinningXMLV17::Class(),this);
   } else {
      R__b.WriteClassBuffer(TUnfoldBinningXMLV17::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TUnfoldBinningXMLV17(void *p) {
      return  p ? new(p) ::TUnfoldBinningXMLV17 : new ::TUnfoldBinningXMLV17;
   }
   static void *newArray_TUnfoldBinningXMLV17(Long_t nElements, void *p) {
      return p ? new(p) ::TUnfoldBinningXMLV17[nElements] : new ::TUnfoldBinningXMLV17[nElements];
   }
   // Wrapper around operator delete
   static void delete_TUnfoldBinningXMLV17(void *p) {
      delete ((::TUnfoldBinningXMLV17*)p);
   }
   static void deleteArray_TUnfoldBinningXMLV17(void *p) {
      delete [] ((::TUnfoldBinningXMLV17*)p);
   }
   static void destruct_TUnfoldBinningXMLV17(void *p) {
      typedef ::TUnfoldBinningXMLV17 current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::TUnfoldBinningXMLV17

namespace {
  void TriggerDictionaryInitialization_TUnfoldV17Dict_Impl() {
    static const char* headers[] = {
0    };
    static const char* includePaths[] = {
"/home/suman/HEP_Package/ROOT616build/include",
"/home/suman/Paradox/Charged_ESV/Working/Unfolding/Ntuple_Fill_RM/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "TUnfoldV17Dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
class __attribute__((annotate(R"ATTRDUMP(Unfolding with support for L-curve analysis)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$TUnfold.h")))  TUnfoldV17;
class __attribute__((annotate(R"ATTRDUMP(iterative Unfolding with scan of SURE)ATTRDUMP"))) TUnfoldIterativeEMV17;
class __attribute__((annotate(R"ATTRDUMP(Unfolding with support for systematic error propagation)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$TUnfoldSys.h")))  TUnfoldSysV17;
class __attribute__((annotate(R"ATTRDUMP(Complex binning schemes for TUnfoldDensity)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$TUnfoldBinning.h")))  TUnfoldBinningV17;
class __attribute__((annotate(R"ATTRDUMP(Unfolding with density regularisation)ATTRDUMP"))) TUnfoldDensityV17;
class __attribute__((annotate(R"ATTRDUMP(Complex binning schemes for TUnfoldDensity)ATTRDUMP"))) TUnfoldBinningXMLV17;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "TUnfoldV17Dict dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
// Author: Stefan Schmitt
// DESY, 19/10/11

//  Version 17.9

//////////////////////////////////////////////////////////////////////////
//                                                                      //
//                                                                      //
//  TUnfold provides functionality to correct data                      //
//   for migration effects.                                             //
//                                                                      //
//  Citation: S.Schmitt, JINST 7 (2012) T10003 [arXiv:1205.6201]        //
//                                                                      //
//  this class implements the iterative EM unfolding method             //
//   (also called D'Agostini Method or "iterative Bayesian method")     //
//  which has been "invented" independently by numnerous authors        //
//  to unfold Poisson-distributed, mutually exclusive bins.             //
//  See e.g.                                                            //
//    Richardson, W.H., Opt. Soc. Amer. A62 (1972), 55                  //
//    Lucy, L.B., Astron. J. 79 (1974), 745.                            //
//    Vardi, Y., Shepp, L.A. and Kaufman, L.,                           //
//         J. Amer. Stat. Assoc. 80 (1985), 8.                          //
//    Multhei, H.N. and Schorr, B., Nucl. Instr. Meth. A257 (1987), 371 //
//    D'Agostini, G.,  Nucl. Instr. Meth. A362 (1995), 487              //
//                                                                      //
//  The novelty with this implementation is that the number of          //
//  iterations can be chosen based on SURE                              //
//        (Stein's unbiased Risk Estimator)                             //
//  See:                                                                //
//    Tibshirani, R.J. and Rosset, S., J. Amer. Stat. Assoc. 114, 526   //
//      [arXiv:1612.09415]                                              //
//                                                                      //
// This method is there for comparison with the Tihkonov unfolding.     //
// The interface is similar to "TUnfoldDensity"                         //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#ifndef ROOT_TUnfoldInterativeEM
#define ROOT_TUnfoldInterativeEM

#include "TUnfold.h"

class TUnfoldBinningV17;

#define TUnfoldIterativeEM TUnfoldIterativeEMV17

class TUnfoldIterativeEMV17 : public TObject {
 public:
   TUnfoldIterativeEMV17(void);
   TUnfoldIterativeEMV17(const TH2 *hist_A, TUnfold::EHistMap histmap,
                         const TUnfoldBinningV17 *outputBins=0,
                         const TUnfoldBinningV17 *inputBins=0);
   virtual ~TUnfoldIterativeEMV17();
   virtual void DoUnfold(Int_t numIterations);
   virtual Int_t SetInput(const TH1 *hist_y, Double_t scaleBias=1.0);
   void SubtractBackground(const TH1 *hist_bgr,const char *name,
                           Double_t scale=1.0);
   void DoUnfold(Int_t nIter,const TH1 *hist_y, Double_t scaleBias=1.0);
   virtual Int_t ScanSURE(Int_t nIterMax,
                          TGraph **SURE=0,
                          TGraph **df_deviance=0);
   TH1 *GetOutput(const char *histogramName,
                  const char *histogramTitle=0,const char *distributionName=0,
		  const char *projectionMode=0,Bool_t useAxisBinning=kTRUE) const;
   TH1 *GetFoldedOutput(const char *histogramName,
                        const char *histogramTitle=0,const char *distributionName=0,
                        const char *projectionMode=0,Bool_t useAxisBinning=kTRUE,
                        Bool_t addBgr=kFALSE) const;
   Double_t GetDeviance(void) const;
   Double_t GetDF(void) const;
   Double_t GetSURE(void) const;
 protected:
   virtual void Reset(void);
   virtual void IterateOnce(void);
   TUnfoldBinningV17 *f_inputBins,*f_outputBins;
   const TUnfoldBinningV17 *f_constInputBins,*f_constOutputBins;
   TMatrixD *fA;
   TVectorD *fEpsilon;
   TVectorD *fX0;
   TVectorD *fY;
   TVectorD *fBgr;
   double fScaleBias;

   Int_t fStep;
   TVectorD *fX;
   TMatrixD *fDXDY;

   ClassDef(TUnfoldIterativeEMV17, TUnfold_CLASS_VERSION) //iterative Unfolding with scan of SURE
};

#endif
// Author: Stefan Schmitt
// DESY, 11/08/11

//  Version 17.9, parallel to changes in TUnfold
//
//  History:
//    Version 17.8, new method GetDXDY()
//    Version 17.7, with bug-fix for curvature regularisation
//    Version 17.6, with updated doxygen comments and bug-fixes in TUnfoldBinning
//    Version 17.5, bug fix in TUnfold also corrects GetEmatrixSysUncorr()
//    Version 17.4, in parallel to changes in TUnfoldBinning
//    Version 17.3, in parallel to changes in TUnfoldBinning
//    Version 17.2, in parallel to changes in TUnfoldBinning
//    Version 17.1, add scan type RhoSquare
//    Version 17.0, support for density regularisation and complex binning schemes

#ifndef ROOT_TUnfoldDensity
#define ROOT_TUnfoldDensity

//////////////////////////////////////////////////////////////////////////
//                                                                      //
//                                                                      //
//  TUnfoldDensity, an extension of the class TUnfoldSys to correct for //
//  migration effects. TUnfoldDensity provides methods to deal with     //
//  multidimensional complex binning schemes and variable bin widths    //
//                                                                      //
//  Citation: S.Schmitt, JINST 7 (2012) T10003 [arXiv:1205.6201]        //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

/*
  This file is part of TUnfold.

  TUnfold is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  TUnfold is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with TUnfold.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "TUnfoldSys.h"
#include "TUnfoldBinning.h"

#define TUnfoldDensity TUnfoldDensityV17

class TUnfoldDensityV17 : public TUnfoldSysV17 {
 protected:
   /// binning scheme for the output (truth level)
   const TUnfoldBinningV17 * fConstOutputBins;
   /// binning scheme for the input (detector level)
   const TUnfoldBinningV17 * fConstInputBins;
   /// pointer to output binning scheme if owned by this class
   TUnfoldBinningV17 *fOwnedOutputBins;
   /// pointer to input binning scheme if owned by this class
   TUnfoldBinningV17 *fOwnedInputBins;
   /// binning scheme for the regularisation conditions
   TUnfoldBinningV17 *fRegularisationConditions;

 public:
   /// choice of regularisation scale factors to cinstruct the matrix L
   enum EDensityMode {
      /// no scale factors, matrix L is similar to unity matrix
      kDensityModeNone=0,
      /// scale factors from multidimensional bin width
      kDensityModeBinWidth=1,
      /// scale factors from user function in TUnfoldBinning
      kDensityModeUser=2,
      /// scale factors from multidimensional bin width and  user function
      kDensityModeBinWidthAndUser=3
   };
 protected:

   virtual TString GetOutputBinName(Int_t iBinX) const; // name a bin

   Double_t GetDensityFactor(EDensityMode densityMode,Int_t iBin) const; // density correction factor for this bin
   void RegularizeDistributionRecursive
     (const TUnfoldBinningV17 *binning,ERegMode regmode,
      EDensityMode densityMode,const char *distribution,
      const char *axisSteering); // regularize the given binning recursively
   void RegularizeOneDistribution
     (const TUnfoldBinningV17 *binning,ERegMode regmode,
      EDensityMode densityMode,const char *axisSteering); // regularize the distribution of one binning node

 public:
   TUnfoldDensityV17(void); // constructor for derived classes, do nothing

   TUnfoldDensityV17(const TH2 *hist_A, EHistMap histmap,
		     ERegMode regmode = kRegModeCurvature,
		     EConstraint constraint=kEConstraintArea,
		     EDensityMode densityMode=kDensityModeBinWidthAndUser,
		     const TUnfoldBinningV17 *outputBins=0,
		     const TUnfoldBinningV17 *inputBins=0,
		     const char *regularisationDistribution=0,
		     const char *regularisationAxisSteering="*[UOB]"); // constructor for using the histogram classes. Default regularisation is on the curvature of the bin-width normalized density, excluding underflow and overflow bins

   virtual ~ TUnfoldDensityV17(void); // delete data members

   void RegularizeDistribution(ERegMode regmode,EDensityMode densityMode,
			       const char *distribution,
			       const char *axisSteering); // regularize distribution(s) of the output binning scheme

   /// scan mode for correlation scan
   enum EScanTauMode {
      /// average global correlation coefficient (from TUnfold::GetRhoI())
      kEScanTauRhoAvg =0,
      /// maximum global correlation coefficient (from TUnfold::GetRhoI())
      kEScanTauRhoMax =1,
       /// average global correlation coefficient (from TUnfoldSys::GetRhoItotal())
      kEScanTauRhoAvgSys =2,
      /// maximum global correlation coefficient (from TUnfoldSys::GetRhoItotal())
      kEScanTauRhoMaxSys =3,
      /// average global correlation coefficient squared (from TUnfold::GetRhoI())
      kEScanTauRhoSquareAvg =4,
      /// average global correlation coefficient squared (from TUnfoldSys::GetRhoItotal())
      kEScanTauRhoSquareAvgSys =5
   };

   virtual Int_t ScanTau(Int_t nPoint,Double_t tauMin,Double_t tauMax,
                         TSpline **scanResult,Int_t mode=kEScanTauRhoAvg,
                         const char *distribution=0,const char *projectionMode=0,TGraph **lCurvePlot=0,TSpline **logTauXPlot=0,TSpline **logTauYPlot=0); // scan some variable (e.g. global correlation) and find a minimum using successive calls to DoUnfold(Double_t) at various tau
   virtual Double_t GetScanVariable(Int_t mode,const char *distribution,const char *projectionMode); // calculate variable for ScanTau()

   TH1 *GetOutput(const char *histogramName,
                  const char *histogramTitle=0,const char *distributionName=0,
		  const char *projectionMode=0,Bool_t useAxisBinning=kTRUE) const;  // get unfolding result
   TH1 *GetBias(const char *histogramName,
                const char *histogramTitle=0,const char *distributionName=0,
		const char *projectionMode=0,Bool_t useAxisBinning=kTRUE) const;      // get bias
   TH1 *GetFoldedOutput(const char *histogramName,
                        const char *histogramTitle=0,
                        const char *distributionName=0,
                        const char *projectionMode=0,Bool_t useAxisBinning=kTRUE,
                        Bool_t addBgr=kFALSE) const; // get unfolding result folded back
   TH1 *GetBackground(const char *histogramName,const char *bgrSource=0,
                      const char *histogramTitle=0,
                      const char *distributionName=0,
                      const char *projectionMode=0,Bool_t useAxisBinning=kTRUE,Int_t includeError=3) const; // get background source
   TH1 *GetInput(const char *histogramName,const char *histogramTitle=0,
                 const char *distributionName=0,
                 const char *projectionMode=0,Bool_t useAxisBinning=kTRUE) const;     // get unfolding input
   TH1 *GetDeltaSysSource(const char *source,
                          const char *histogramName,
                          const char *histogramTitle=0,
                          const char *distributionName=0,
			  const char *projectionMode=0,Bool_t useAxisBinning=kTRUE); // get systematic shifts from one systematic source
   TH1 *GetDeltaSysBackgroundScale(const char *bgrSource,
                                   const char *histogramName,
                                   const char *histogramTitle=0,
                                   const char *distributionName=0,
				   const char *projectionMode=0,Bool_t useAxisBinning=kTRUE); // get correlated uncertainty induced by the scale uncertainty of a background source
   TH1 *GetDeltaSysTau(const char *histogramName,
                       const char *histogramTitle=0,
                       const char *distributionName=0,
		       const char *projectionMode=0,Bool_t useAxisBinning=kTRUE); // get correlated uncertainty from varying tau
   TH2 *GetEmatrixSysUncorr(const char *histogramName,
			    const char *histogramTitle=0,
			    const char *distributionName=0,
			    const char *projectionMode=0,Bool_t useAxisBinning=kTRUE); // get error matrix contribution from uncorrelated errors on the matrix A
   TH2 *GetEmatrixSysBackgroundUncorr(const char *bgrSource,
				      const char *histogramName,
				      const char *histogramTitle=0,
				      const char *distributionName=0,
				      const char *projectionMode=0,Bool_t useAxisBinning=kTRUE); // get error matrix from uncorrelated error of one background source
   TH2 *GetEmatrixInput(const char *histogramName,
                        const char *histogramTitle=0,
			const char *distributionName=0,
			const char *projectionMode=0,Bool_t useAxisBinning=kTRUE); // get error contribution from input vector
   TH2 *GetEmatrixTotal(const char *histogramName,
			const char *histogramTitle=0,
			const char *distributionName=0,
			const char *projectionMode=0,Bool_t useAxisBinning=kTRUE); // get total error including systematic,statistical,background,tau errors
   TH1 *GetRhoIstatbgr(const char *histogramName,const char *histogramTitle=0,
                     const char *distributionName=0,
		     const char *projectionMode=0,Bool_t useAxisBinning=kTRUE,
                     TH2 **ematInv=0);      // get global correlation coefficients, stat+bgr errors only (from TUnfold)
   TH1 *GetRhoItotal(const char *histogramName,const char *histogramTitle=0,
                     const char *distributionName=0,
		     const char *projectionMode=0,Bool_t useAxisBinning=kTRUE,
                     TH2 **ematInv=0);      // get global correlation coefficients, including systematic errors (from TUnfoldSys)
   TH2 *GetRhoIJtotal(const char *histogramName,
		      const char *histogramTitle=0,
		      const char *distributionName=0,
		      const char *projectionMode=0,Bool_t useAxisBinning=kTRUE);     // get correlation coefficients
   TH2 *GetL(const char *histogramName,
             const char *histogramTitle=0,
             Bool_t useAxisBinning=kTRUE); // get regularisation matrix
   TH1 *GetLxMinusBias(const char *histogramName,const char *histogramTitle=0); // get vector L(x-bias) of regularisation conditions 

   TH2 *GetProbabilityMatrix(const char *histogramName,
                           const char *histogramTitle=0,Bool_t useAxisBinning=kTRUE) const; // get matrix of probabilities
   TH2 *GetDXDY(const char *histogramName,
                          const char *histogramTitle=0,Bool_t useAxisBinning=kTRUE) const; // get matrix DX/DY

   const TUnfoldBinningV17 *GetInputBinning(const char *distributionName=0) const; // find binning scheme for input bins
   const TUnfoldBinningV17 *GetOutputBinning(const char *distributionName=0) const; // find binning scheme for output bins
   /// return binning scheme for regularisation conditions (matrix L)
TUnfoldBinningV17 *GetLBinning(void) const { return fRegularisationConditions; }
   ClassDef(TUnfoldDensityV17, TUnfold_CLASS_VERSION) //Unfolding with density regularisation
};

#endif
// Author: Stefan Schmitt
// DESY, 10/08/11

//  Version 17.9, parallel to changes in TUnfold
//
//  History:
//    Version 17.8, relaxed DTD definition
//    Version 17.7, in parallel to changes in TUnfold
//    Version 17.6, with updated doxygen comments
//    Version 17.5, in parallel to changes in TUnfold
//    Version 17.4, in parallel to changes in TUnfoldBinning
//    Version 17.3, support for repeated bins with the same width
//    Version 17.2, XML interface for class TUnfoldBinning

#ifndef ROOT_TUnfoldBinningXML
#define ROOT_TUnfoldBinningXML


//////////////////////////////////////////////////////////////////////////
//                                                                      //
//                                                                      //
//  TUnfoldBinningXML, an auxillary class to read and write             //
//  complex binning schemes in XML                                      //
//                                                                      //
//  Citation: S.Schmitt, JINST 7 (2012) T10003 [arXiv:1205.6201]        //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

/*
  This file is part of TUnfold.

  TUnfold is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  TUnfold is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with TUnfold.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "TUnfoldBinning.h"
#include <iostream>
#include <TNamed.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TXMLNode.h>
#include <TXMLDocument.h>
#include <iostream>

class TXMLNode;
class TXMLDocument;

#define TUnfoldBinningXML TUnfoldBinningXMLV17

class TUnfoldBinningXMLV17 : public TUnfoldBinningV17 {
 public:
   /********************** XML interface to read binning schemes *************/
static TUnfoldBinningXMLV17 *ImportXML(const TXMLDocument *document,const char *name); // import binning scheme
   static Int_t ExportXML(const TUnfoldBinningV17 &binning,std::ostream &out,Bool_t writeHeader,Bool_t writeFooter,Int_t indent=0); // append binning scheme to file
   Int_t ExportXML(const char *fileName) const; // export this binning scheme
   static void WriteDTD(const char *fileName="tunfoldbinning.dtd"); // write DTD file
   static void WriteDTD(std::ostream &out); // write DTD to stream

   /// construct a new binning scheme, for use with the root streamer
   TUnfoldBinningXMLV17 (const char *name=0,Int_t nBins=0,const char *binNames=0) 
      : TUnfoldBinningV17 (name,nBins,binNames) { }
 protected:
   static TUnfoldBinningXMLV17 *ImportXMLNode(TXMLNode *node); // import the given node as binning scheme
   void AddAxisXML(TXMLNode *node); // import axis information
protected:

   ClassDef(TUnfoldBinningXMLV17, TUnfold_CLASS_VERSION) //Complex binning schemes for TUnfoldDensity
};

#endif

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"TUnfoldBinningV17", payloadCode, "@",
"TUnfoldBinningXMLV17", payloadCode, "@",
"TUnfoldDensityV17", payloadCode, "@",
"TUnfoldIterativeEMV17", payloadCode, "@",
"TUnfoldSysV17", payloadCode, "@",
"TUnfoldV17", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("TUnfoldV17Dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_TUnfoldV17Dict_Impl, {}, classesHeaders, /*has no C++ module*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_TUnfoldV17Dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_TUnfoldV17Dict() {
  TriggerDictionaryInitialization_TUnfoldV17Dict_Impl();
}
