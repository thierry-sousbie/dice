#ifndef __HLP_TYPE_LIST_HXX__
#define __HLP_TYPE_LIST_HXX__

#include "../../internal/namespace.header"

namespace hlp
{

	class EmptyList
	{
	};

	template <class F1 = EmptyList, class F2 = EmptyList, class F3 = EmptyList,
			  class F4 = EmptyList, class F5 = EmptyList, class F6 = EmptyList,
			  class F7 = EmptyList, class F8 = EmptyList, class F9 = EmptyList,
			  class F10 = EmptyList, class F11 = EmptyList, class F12 = EmptyList,
			  class F13 = EmptyList, class F14 = EmptyList, class F15 = EmptyList,
			  class F16 = EmptyList, class F17 = EmptyList, class F18 = EmptyList,
			  class F19 = EmptyList>
	class TypeListT;

	template <>
	class TypeListT<EmptyList, EmptyList, EmptyList,
					EmptyList, EmptyList, EmptyList,
					EmptyList, EmptyList, EmptyList,
					EmptyList, EmptyList, EmptyList,
					EmptyList, EmptyList, EmptyList,
					EmptyList, EmptyList, EmptyList,
					EmptyList>
	{
	public:
		typedef TypeListT<EmptyList, EmptyList, EmptyList,
						  EmptyList, EmptyList, EmptyList,
						  EmptyList, EmptyList, EmptyList,
						  EmptyList, EmptyList, EmptyList,
						  EmptyList, EmptyList, EmptyList,
						  EmptyList, EmptyList, EmptyList,
						  EmptyList>
			MyType;
		typedef MyType Next;
		typedef EmptyList Type;

		static const int INDEX = -1;
		static const int SIZE = INDEX + 1;
	};

	template <class F1, class F2, class F3,
			  class F4, class F5, class F6,
			  class F7, class F8, class F9,
			  class F10, class F11, class F12,
			  class F13, class F14, class F15,
			  class F16, class F17, class F18,
			  class F19>
	class TypeListT : public TypeListT<F2, F3, F4, F5, F6, F7, F8, F9, F10,
									   F11, F12, F13, F14, F15, F16, F17, F18, F19>
	{
	public:
		typedef TypeListT<F1, F2, F3, F4, F5, F6, F7, F8, F9, F10,
						  F11, F12, F13, F14, F15, F16, F17, F18, F19>
			MyType;
		typedef TypeListT<F2, F3, F4, F5, F6, F7, F8, F9, F10,
						  F11, F12, F13, F14, F15, F16, F17, F18, F19>
			Next;
		typedef F1 Type;

		static const int INDEX = (Next::INDEX + 1);
		static const int SIZE = INDEX + 1;
	};

}

#include "../../internal/namespace.footer"
#endif
