#include <iostream>
#include <string>
#include <cstring>
#include <stdexcept>
#include <vector>
#include <memory>
#include <assert.h>
#include <fstream>
#include <streambuf>


#define MAX_VALUE_IN_BITS(X) ((1 << X) - 1)

#define LEN_BITS 8
#define MAX_WORD (MAX_VALUE_IN_BITS(LEN_BITS) - 1)
#define USE_STATIC_LEN false
#define WORD_TRESHOLD 2
#define USE_WORD_TRESHOLD false
#define USE_DYNAMIC_PRIORITY false
#define USE_SYMBOL_PAIRS false
#define USE_TABLE_HEADER false
#define USE_LONG_TABLE false
#define USE_UNICODE false
#define USE_CYRILLIC_PAIRS true
#define USE_RLE false
#define RLE_BITS 7
#define RLE_MAX_LEN (MAX_VALUE_IN_BITS(RLE_BITS))
#define USE_RLE_DYNAMIC_LENGTH false




using std::string;

class CBitOStream {
public:
    CBitOStream();
    void WriteBit( bool bit );
    void WriteByte( unsigned char byte );
    // Закрывает поток - перебрасывает буфер в результирующую строку, дописывает число бит в последнем байте.
    void Close();
    
    // Получить результат можно только после Close.
    string GetResult() const;
    
private:
    string result;
    unsigned char buffer;
    char bitsInBuffer;
    bool isClosed;
};

CBitOStream::CBitOStream() :
buffer( 0 ),
bitsInBuffer( 0 ),
isClosed( false )
{
}

void CBitOStream::WriteBit( bool bit )
{
    assert( !isClosed );
    
    // Устанавливаем или сбрасываем нужный бит. Первыми битами считаем старшие.
    if( bit ) {
        buffer |= 1 << ( 7 - bitsInBuffer );
    } else {
        buffer &= ~( 1 << ( 7 - bitsInBuffer ) );
    }
    ++bitsInBuffer;
    if( bitsInBuffer == 8 ) {
        // Буфер закончился. Сбрасываем в результирующую строку.
        result.push_back( buffer );
        buffer = 0;
        bitsInBuffer = 0;
    }
}

void CBitOStream::WriteByte( unsigned char byte )
{
    assert( !isClosed );
    
    buffer |= byte >> ( bitsInBuffer ); // Пользуемся тем, что оставшиеся с buffer байты нулевые.
    unsigned char nextBuffer = byte << ( 8 - bitsInBuffer );
    
    // Сбрасываем буфер в результирующую строку.
    result.push_back( buffer );
    buffer = nextBuffer;
}

void CBitOStream::Close()
{
    assert( !isClosed );
    
    if( bitsInBuffer == 0 ) {
        // Дописываем, что в последнем байте 8 значащих бит.
        result.push_back( 8 );
    } else {
        // Скидываем недозаполненный буфер и количество значащих бит.
        result.push_back( buffer );
        result.push_back( bitsInBuffer );
    }
    isClosed = true;
}

string CBitOStream::GetResult() const
{
    assert( isClosed );
    return result;
}

// ----------------------------------------------------------------------------------------------------------

class CBitIStream {
public:
    explicit CBitIStream( const string& _source );
    
    // Чтение неоконченного потока.
    bool ReadBit();
    unsigned char ReadByte();
    
    bool IsFinished() const;
    
private:
    const string source;
    int bitsCount; // Общее число бит.
    int pos; // Считанное число бит.
};

CBitIStream::CBitIStream( const string& _source ) :
source( _source ),
pos( 0 )
{
    assert( !source.empty() );
    
    // Общее число бит.
    bitsCount = 8 * ( (int)source.length() - 2 ) + source[source.length() - 1];
}


bool CBitIStream::ReadBit()
{
    unsigned int mask = 1 << (7 - (pos % 8));
    bool value = (((unsigned char)(source[pos / 8])) & mask) != 0;
    pos++;
    return value;
}


unsigned char CBitIStream::ReadByte()
{
    unsigned char value = 0;
    
    for (int i = 0; i < 8; i++) {
        unsigned char bit = (ReadBit()) ? (1) : (0);
        bit = bit << (7 - i);
        value |= bit;
    }
    
    return value;
}


bool CBitIStream::IsFinished() const
{
    return pos >= bitsCount;
}






#define INIT_SIZE 8
#define INC_FACTOR 0.75f
#define DEC_FACTOR 0.25f



typedef unsigned int Hash;
template<typename T> using HashFunction = Hash(*)(const T& data, const int len, bool arg);
Hash gornerCharHash(char* data, int len, int dataLen, bool arg);
template<typename T> using Deallocator = void(*)(T& data);
template<typename T> using Allocator = T(*)(T& copyFrom);
template<typename T> using Comparator = bool(*)(T& a, T& b);
template<typename T> using LessComparator = bool(*)(T& a, T& b);






template<typename T> Hash GornerHash(const T& data, const int len, bool arg)
{
    return gornerCharHash((char*) &data, sizeof(T), len, arg);
}


template<> Hash GornerHash<char*>(char* const& data, const int len, bool arg)
{
    return gornerCharHash((char*) data, sizeof(char) * (int) strlen(data), len, arg);
}


template<typename T> void DefaultDeallocator(T& data)
{
}


template<typename T> T DefaultAllocator(T& copyFrom)
{
    return copyFrom;
}


template<typename T> bool DefaultComparator(T& a, T& b)
{
    return a == b;
}


template<typename T> bool DefaultLessComparator(T& a, T& b)
{
    return a < b;
}


template<> void DefaultDeallocator<char*>(char*& data)
{
    delete[] ((char*) data);
}


template<> char* DefaultAllocator<char*>(char*& copyFrom)
{
    if (copyFrom == NULL) {
        char* res = new char[1];
        res[0] = 0;
        return res;
    }
    
    int len = (int) strlen(copyFrom);
    char* res = new char[len + 1];
    strcpy(res, copyFrom);
    res[len] = 0;
    return res;
}


template<> bool DefaultComparator<char*>(char*& a, char*& b)
{
    if ((a == NULL) || (b == NULL)) return false;
    return strcmp((char*) a, (char*) b) == 0;
}



Hash gornerCharHash(char* data, int len, int dataLen, bool arg)
{
    Hash a = 1073676287;
    if (arg) a = 2764443;
    if (len <= 0) return 0;
    
    Hash res = 0;
    
    for (int i = 0; i < len; i++) {
        res = ((res * a) + data[i]) % dataLen;
    }
    
    return res;
}




template<
typename K, typename V,
HashFunction<K> HashFunction = GornerHash<K>,
Allocator<K> KeyAllocator = DefaultAllocator<K>,
Allocator<V> ValueAllocator = DefaultAllocator<V>,
Deallocator<K> KeyDeallocator = DefaultDeallocator<K>,
Deallocator<V> ValueDeallocator = DefaultDeallocator<V>,
Comparator<K> KeyComparator = DefaultComparator<K>
> class HashMap
{
public:
    struct Pair {
        K key;
        V value;
        bool empty;
        bool deleted = false;
        
        const bool isEmpty() { return empty || deleted; }
        
        Pair(K& key, V& value, bool empty = true)
        { this->key = key; this->value = value; this->empty = empty; }
        
        Pair() { empty = true; }
    };
    
    struct Iter {
        Pair* data;
        int size;
        int pos;
        
        Iter() {}
        Iter(const Iter& iter): Iter(iter.data, iter.size, iter.pos) {}
        Iter(Pair* data, int size, int pos)
        : data(data), size(size)
        { this->pos = pos; if (pos <= 0) { this->pos = -1; this->operator++();} }
        
        void operator++();
        const bool operator!=(const Iter& iter)
        { return this->pos != iter.pos; }
    };
    
    struct PairIter: public Iter {
        PairIter() {}
        PairIter(const PairIter& iter): PairIter(iter.data, iter.size, iter.pos) {}
        PairIter(Pair* data, int size, int pos)
        : Iter(data, size, pos) {}
        
        Pair& operator*()
        {
            if (this->pos < this->size) return this->data[this->pos];
            throw std::runtime_error("Dereferencing end() iterator");
        }
    };
    
    struct KeyIter: public Iter {
        KeyIter() {}
        KeyIter(const KeyIter& iter): KeyIter(iter.data, iter.size, iter.pos) {}
        KeyIter(Pair* data, int size, int pos)
        : Iter(data, size, pos) {}
        
        K& operator*()
        {
            if (this->pos < this->size) return this->data[this->pos].key;
            throw std::runtime_error("Dereferencing end() iterator");
        }
    };
    
    struct ValueIter: public Iter {
        ValueIter() {}
        ValueIter(const Iter& iter): ValueIter(iter.data, iter.size, iter.pos) {}
        ValueIter(Pair* data, int size, int pos)
        : Iter(data, size, pos) {}
        
        V& operator*()
        {
            if (this->pos < this->size) return this->data[this->pos].value;
            throw std::runtime_error("Dereferencing end() iterator");
        }
    };
    
private:
    template<typename I = PairIter> struct Iters {
        I _begin;
        I _end;
        
        const I& begin() { return _begin; }
        const I& end() { return _end; }
        
        Iters(Pair* data, int size)
        {
            _begin = I(data, size, 0);
            _end = I(data, size, size);
        }
    };
    
private:
    Pair* data = NULL;
    int size = 0;
    int used = 0;
    int minSize = 0;
    
protected:
    const Hash getHash(const K& key, const int i);
    void resize(int size);
    void clean();
    float loadFactor() { return (float) used / (float) size; }
    
public:
    HashMap(int init = INIT_SIZE) { minSize = init; resize(init); }
    ~HashMap() { clean(); }
    
    V get(K key);
    void set(K key, V value, bool rewrite = true);
    const bool has(K key);
    void remove(K key);
    
    V get(K key, V def)
    {
        try {
            return get(key);
        } catch(...) {
            return def;
        }
    }
    
    void put(K key, V value) { set(key, value, false); }
    
    void clear() {
        clean();
        resize(minSize);
    }
    
    Iters<PairIter> pairs() { return Iters<PairIter>(data, size); }
    Iters<KeyIter> keys() { return Iters<KeyIter>(data, size); }
    Iters<ValueIter> values() { return Iters<ValueIter>(data, size); }
    
    PairIter begin() { return pairs().begin(); }
    PairIter end() { return pairs().end(); }
    
    const int count() { return used; }
    
    void debug()
    {
        for (int i = 0; i < size; i++) {
            if (data[i].empty) {
                std::cout << i << ", NULL" << std::endl;
            } else {
                std::cout << i << ", " << data[i].key;
                std::cout << " = " << data[i].value << std::endl;
            }
        }
    }
};




template<
typename K, typename V,
HashFunction<K> HashFunction,
Allocator<K> KeyAllocator,
Allocator<V> ValueAllocator,
Deallocator<K> KeyDeallocator,
Deallocator<V> ValueDeallocator,
Comparator<K> KeyComparator
> void HashMap<K, V, HashFunction, KeyAllocator, ValueAllocator, KeyDeallocator, ValueDeallocator, KeyComparator>
::set(K key, V value, bool rewrite)
{
    Hash h1 = getHash(key, 0);
    Hash h2 = getHash(key, 1);
    Hash h = 0;
    int i = -1;
    int firstDel = -1;
    
    for (i = 0; i < size; i++) {
        h = (h1 + i * h2) % size;
        
        if ((firstDel == -1) && (data[h].deleted)) firstDel = h;
        if (data[h].empty) break;
        
        if (!data[h].deleted && KeyComparator(key, data[h].key)) {
            if (!rewrite) throw std::runtime_error("Duplicate key found");
            used--;
            break;
        }
    }
    
    if (i == size) {
        h = firstDel;
        KeyDeallocator(data[h].key);
        ValueDeallocator(data[h].value);
    }
    
    if (i == -1) throw std::runtime_error("HashMap internal error");
    
    K keyAlloc= KeyAllocator(key);
    V valueAlloc = ValueAllocator(value);
    data[h] = Pair(keyAlloc, valueAlloc, false);
    used++;
    
    if (loadFactor() >= INC_FACTOR) {
        resize(size * 2);
    }
}



template<
typename K, typename V,
HashFunction<K> HashFunction,
Allocator<K> KeyAllocator,
Allocator<V> ValueAllocator,
Deallocator<K> KeyDeallocator,
Deallocator<V> ValueDeallocator,
Comparator<K> KeyComparator
> void HashMap<K, V, HashFunction, KeyAllocator, ValueAllocator, KeyDeallocator, ValueDeallocator, KeyComparator>
::clean() {
    
    if (data == NULL) return;
    
    for (Pair& pair: this->pairs()) {
        if (!pair.empty) {
            KeyDeallocator(pair.key);
            ValueDeallocator(pair.value);
        }
    }
    
    delete[] data;
    data = NULL;
    size = 0;
    used = 0;
}




template<
typename K, typename V,
HashFunction<K> HashFunction,
Allocator<K> KeyAllocator,
Allocator<V> ValueAllocator,
Deallocator<K> KeyDeallocator,
Deallocator<V> ValueDeallocator,
Comparator<K> KeyComparator
> void HashMap<K, V, HashFunction, KeyAllocator, ValueAllocator, KeyDeallocator, ValueDeallocator, KeyComparator>
::resize(int size) {
    
    if (data == NULL) {
        data = new Pair[size];
        this->size = size;
        return;
    }
    
    if (size < minSize) size = minSize;
    HashMap newMap(size);
    newMap.minSize = minSize;
    
    for (Pair& pair: this->pairs()) {
        newMap.put(pair.key, pair.value);
    }
    
    clean();
    this->data = newMap.data;
    this->used = newMap.used;
    this->size = newMap.size;
    newMap.data = NULL;
}




template<
typename K, typename V,
HashFunction<K> HashFunction,
Allocator<K> KeyAllocator,
Allocator<V> ValueAllocator,
Deallocator<K> KeyDeallocator,
Deallocator<V> ValueDeallocator,
Comparator<K> KeyComparator
> V HashMap<K, V, HashFunction, KeyAllocator, ValueAllocator, KeyDeallocator, ValueDeallocator, KeyComparator>
::get(K key)
{
    Hash h1 = getHash(key, 0);
    Hash h2 = getHash(key, 1);
    Hash h = 0;
    int i;
    
    for (i = 0; i < size; i++) {
        h = (h1 + i * h2) % size;
        
        if (data[h].empty) break;
        if (!data[h].deleted && KeyComparator(key, data[h].key))
            return data[h].value;
    }
    
    throw std::runtime_error("Key not found");
}




template<
typename K, typename V,
HashFunction<K> HashFunction,
Allocator<K> KeyAllocator,
Allocator<V> ValueAllocator,
Deallocator<K> KeyDeallocator,
Deallocator<V> ValueDeallocator,
Comparator<K> KeyComparator
> const bool HashMap<K, V, HashFunction, KeyAllocator, ValueAllocator, KeyDeallocator, ValueDeallocator, KeyComparator>
::has(K key)
{
    
    Hash h1 = getHash(key, 0);
    Hash h2 = getHash(key, 1);
    Hash h = 0;
    int i;
    
    for (i = 0; i < size; i++) {
        h = (h1 + i * h2) % size;
        
        if (data[h].empty) return false;
        if (!data[h].deleted && KeyComparator(key, data[h].key))
            return true;
    }
    
    return false;
}




template<
typename K, typename V,
HashFunction<K> HashFunction,
Allocator<K> KeyAllocator,
Allocator<V> ValueAllocator,
Deallocator<K> KeyDeallocator,
Deallocator<V> ValueDeallocator,
Comparator<K> KeyComparator
> void HashMap<K, V, HashFunction, KeyAllocator, ValueAllocator, KeyDeallocator, ValueDeallocator, KeyComparator>
::remove(K key)
{
    
    Hash h1 = getHash(key, 0);
    Hash h2 = getHash(key, 1);
    Hash h = 0;
    int i;
    
    for (i = 0; i < size; i++) {
        h = (h1 + i * h2) % size;
        
        if (data[h].empty) throw std::runtime_error("Key not found");
        
        if (!data[h].deleted && KeyComparator(key, data[h].key)) {
            data[h].deleted = true;
            used--;
            
            if (loadFactor() <= DEC_FACTOR) {
                resize(used / 2);
            }
            
            return;
        }
    }
    
    throw std::runtime_error("Key not found");
}




template<
typename K, typename V,
HashFunction<K> HashFunction,
Allocator<K> KeyAllocator,
Allocator<V> ValueAllocator,
Deallocator<K> KeyDeallocator,
Deallocator<V> ValueDeallocator,
Comparator<K> KeyComparator
> void HashMap<K, V, HashFunction, KeyAllocator, ValueAllocator, KeyDeallocator, ValueDeallocator, KeyComparator>
::Iter::operator++()
{
    
    for (int i = pos + 1; i < size; i++) {
        if ((i < size) && (!data[i].isEmpty())) {
            pos = i;
            return;
        }
    }
    
    pos = size;
}




template<
typename K, typename V,
HashFunction<K> HashFunction,
Allocator<K> KeyAllocator,
Allocator<V> ValueAllocator,
Deallocator<K> KeyDeallocator,
Deallocator<V> ValueDeallocator,
Comparator<K> KeyComparator
> const Hash HashMap<K, V, HashFunction, KeyAllocator, ValueAllocator, KeyDeallocator, ValueDeallocator, KeyComparator>
::getHash(const K &key, const int i)
{
    if (i <= 0) return HashFunction(key, size, false) % size;
    
    Hash res = HashFunction(key, size, true);
    res = (res * 2 + 1) % size;
    return res;
}










template<
typename T,
Comparator<T> Comparator = DefaultComparator<T>,
LessComparator<T> LessComparator = DefaultLessComparator<T>
> class Heap
{
private:
    std::vector<T> buff;
    
protected:
    void shiftUp(int index);
    void shiftDown(int index);
    
public:
    Heap() {}
    
    void push(T el);
    T pop();
    
    T get(int index) { return buff[index]; }
    const size_t size() { return buff.size(); }
    
    void clear() { buff.clear(); }
    
    void debug() {
        for (int i = 0; i < (int) buff.size(); i++) {
            std::cout << i << ", " << buff[i] << std::endl;
        }
    }
};



template<
typename T,
Comparator<T> Comparator,
LessComparator<T> LessComparator
> void Heap<T, Comparator, LessComparator>::shiftUp(int index)
{
    
    while (index > 0) {
        int parent = (index - 1) / 2;
        if (!LessComparator(buff[parent], buff[index])) return;
        std::swap(buff[index], buff[parent]);
        index = parent;
        if (index <= 0) return;
    }
}


template<
typename T,
Comparator<T> Comparator,
LessComparator<T> LessComparator
> void Heap<T, Comparator, LessComparator>::shiftDown(int index)
{
    
    bool running = true;
    
    while (running) {
        int left = 2 * index + 1;
        int right = left + 1;
        int largest = left;
        running = false;
        
        if ((left < buff.size()) && (right < buff.size())) {
            if (LessComparator(buff[largest], buff[right])) {
                largest = right;
            }
        }
        
        if (largest < size()) {
            if (LessComparator(buff[index], buff[largest])) {
                std::swap(buff[index], buff[largest]);
                index = largest;
                running = true;
            }
        }
    }
}


template<
typename T,
Comparator<T> Comparator,
LessComparator<T> LessComparator
> void Heap<T, Comparator, LessComparator>::push(T el)
{
    
    buff.push_back(el);
    shiftUp((int) size() - 1);
}


template<
typename T,
Comparator<T> Comparator,
LessComparator<T> LessComparator
> T Heap<T, Comparator, LessComparator>::pop()
{
    
    if (size() == 0)
        throw std::length_error("Attempt to pop from an empty heap");
    
    T res = buff[0];
    buff[0] = buff[size() - 1];
    buff.resize(size() - 1);
    
    if (size() != 0) shiftDown(0);
    return res;
}








typedef unsigned char byte;

struct Node;
struct Bytes;
std::ostream& operator<<(std::ostream& os, const Node& obj);
std::ostream& operator<<(std::ostream& os, const Bytes& obj);


struct Node
{
    byte code = 0;
    bool hasCode = false;
    int count = 0;
    
    Node* left = NULL;
    Node* right = NULL;
    Node* parent = NULL;
    
    std::string key;
    
    void clear()
    { if (left != NULL) delete left; if (right != NULL) delete right; }
    
    Node(std::string key = ""):key(key) {}
    Node(const Node& node)
    {
        code = node.code; hasCode = node.hasCode; count = node.count;
        left = node.left; right = node.right; key = node.key;
        parent = node.parent;
    }
    
    ~Node() {}
    
    void debug(int tabs = 0) {
        if (left != NULL) left->debug(tabs - 1);
        std::cout << tabs << ": " << *this << std::endl;
        if (right != NULL) right->debug(tabs + 1);
    }
};


struct Bytes
{
    std::vector<byte> data;
    
    Bytes() {}
    Bytes(const Bytes& bytes) { data = bytes.data; }
    Bytes(std::vector<byte>& bytes) { data = bytes; }
    ~Bytes() {}
};


std::ostream& operator<<(std::ostream& os, const Node& obj)
{
    os << "{ key = " << obj.key;
    
    if (obj.hasCode) {
        os << ", code = " << (int) obj.code;
    }
    
    os << ", count = " << obj.count << " }";
    return os;
}


std::ostream& operator<<(std::ostream& os, const Bytes& obj)
{
    for (const byte& b: obj.data) os << (int)b;
    return os;
}


template<> bool DefaultComparator(Node& a, Node& b)
{
    return (a.count == b.count) && (a.code == b.code);
}


template<> bool DefaultLessComparator(Node& a, Node& b)
{
    if (USE_DYNAMIC_PRIORITY) {
        int countDelta = a.count - b.count;
        if (countDelta < 0) countDelta *= -1;
        
        int keyDelta = (int) a.key.length() - (int) b.key.length();
        if (keyDelta < 0) keyDelta *= -1;
        
        if ((countDelta < 10) && (keyDelta > 3)) {
            return a.key.length() < b.key.length();
        }
        
        return a.count > b.count;
    }
    
    return a.count > b.count;  // reverse
}


template<> Hash GornerHash<std::string>(std::string const& data, const int len, bool arg)
{
    char* str = (char*) data.c_str();
    return gornerCharHash(str, (int) data.length(), len, arg);
}



class Compressor
{
private:
    std::string src;
    Heap<Node> heap;
    HashMap<std::string, Node> dict;
    std::shared_ptr<Node> tree;
    HashMap<std::string, Bytes> map;
    std::string header;
protected:
    void buildStat();
    void buildHeap();
    void buildTree();
    void buildMap(std::vector<byte> path = std::vector<byte>(),
                  Node* node = NULL);
    void buildHeader(Node* node = NULL);
    
    void writeHeader(CBitOStream& os);
    void writeBits(CBitOStream& os, const Bytes& bytes);
    void writeLen(CBitOStream& os, const int len);
    void writeByteLen(CBitOStream& os, const int len);
    
    void writeByteLenExplicit(CBitOStream &os, const int len, int bits);
    int readByteLenExplicit(CBitIStream &is, int bits);
    
    void writeBitsRLE(CBitOStream &os, const Bytes &bytes);
    int readBitsRLE(CBitIStream &is);
    
    void readHeader(CBitIStream& is);
    Bytes readBits(CBitIStream& is, const int len);
    int readLen(CBitIStream& os);
    int readByteLen(CBitIStream &is);
    
    void parseHeader();
    
    bool isMulti() { return false; /* src.length() > 10;*/ }
    bool isLetter(byte ch);
    bool bitsCompare(const Bytes &a, const Bytes &b);
public:
    Compressor(std::string data): src(data) {}
    ~Compressor() { if (tree.get() != NULL) tree->clear(); }
    
    std::string compressed();
    std::string decompressed();
};





bool Compressor::isLetter(byte ch)
{
    static bool lastUnicode = false;
    
    if (USE_UNICODE) {
        if (ch >= 128) {
            lastUnicode = true;
            return true;
        }
        
        if (lastUnicode) {
            lastUnicode = false;
            return false;
        }
    }
    
    return (((ch >= 48) && (ch <= 57)) ||
            ((ch >= 65) && (ch <= 90)) ||
            ((ch >= 97) && (ch <= 122))) || (ch == '#') || (ch == '_'); /* ||
                                                                         (ch == '+') || (ch == '-');*/
    
    
    /*||
     (ch == '(') || (ch == ')') || (ch == '{') || (ch == '}') ||
     (ch == '[') || (ch == ']') || (ch == '+') || (ch == '-');*/
}


void Compressor::buildStat()
{
    HashMap<std::string, bool> symbolPairs;
    
    if (USE_SYMBOL_PAIRS) {
        symbolPairs.put("::", true);
        symbolPairs.put("->", true);
        symbolPairs.put("++", true);
        symbolPairs.put("--", true);
    }
    
    if (USE_CYRILLIC_PAIRS) {
        std::string small = "йцукенгшщзхъфывапролджэячсмитьбю";
        std::string caps = "ЙЦУКЕНГШЩЗХЪФЫВАПРОЛДЖЭЯЧСМИТЬБЮ";
        
        for (int i = 0; i < (int) small.length() / 2; i++) {
            std::string pair;
            
            pair.push_back(small[2 * i]);
            pair.push_back(small[2 * i + 1]);
            symbolPairs.put(pair, true);
            
            pair.clear();
            pair.push_back(caps[2 * i]);
            pair.push_back(caps[2 * i + 1]);
            symbolPairs.put(pair, true);
        }
    }
    
    std::string word;
    
    for (int i = 0; i < (int) src.length(); i++) {
        if (isLetter(src[i])) {
            word.push_back(src[i]);
        }
        
        if ((!isLetter(src[i])) || (i == (int) src.length() - 1) ||
            (word.length() > MAX_WORD)) {
            if (word.length() > 0) {
                Node node = dict.get(word, Node(word));
                node.count++;
                dict.set(word, node);
                word.clear();
            }
            
            if (!isLetter(src[i])) {
                if (i != (int) src.length() - 1) {
                    std::string sym;
                    sym.push_back(src[i]);
                    sym.push_back(src[i + 1]);
                    
                    if (symbolPairs.has(sym)) {
                        i++;
                        continue;
                    }
                }
                
                std::string sym;
                sym.push_back(src[i]);
                
                Node node = dict.get(sym, Node(sym));
                node.count++;
                dict.set(sym, node);
            }
        }
    }
    
    if (USE_SYMBOL_PAIRS || USE_CYRILLIC_PAIRS) {
        for (auto& pair: symbolPairs.keys()) {
            dict.set(pair, Node(pair));
        }
        
        word.clear();
        
        for (int i = 0; i < (int) src.length(); i++) {
            word.push_back(src[i]);
            if (word.length() > 2) word.erase(word.begin());
            
            if (word.length() == 2) {
                if (symbolPairs.has(word)) {
                    Node node = dict.get(word);
                    node.count++;
                    dict.set(word, node);
                    //std::cout << "found " << word << std::endl;
                    word.clear();
                } else {
                    //std::cout << "did_not_find " << word << std::endl;
                }
            }
        }
    }
    
    if (USE_WORD_TRESHOLD) {
        std::vector<std::string> toDelete;
        HashMap<byte, bool> chars;
        
        for (auto& pair: dict) {
            if ((pair.key.length() > 1) &&
                ((pair.value.count < WORD_TRESHOLD))) {
                toDelete.push_back(pair.key);
                
                for (auto& ch: pair.key) {
                    chars.set(ch, true);
                }
            }
        }
        
        for (std::string& key: toDelete) {
            dict.remove(key);
        }
        
        for (int i = 0; i < (int) src.length(); i++) {
            if (chars.has(src[i])) {
                std::string sym;
                sym.push_back(src[i]);
                
                Node node = dict.get(sym, Node(sym));
                node.count++;
                dict.set(sym, node);
            }
        }
    }
}


void Compressor::buildHeap()
{
    if (USE_TABLE_HEADER) {
        Heap<Node> nodes;
        
        for (Node& node: dict.values()) {
            if (node.count > 0) nodes.push(node);
        }
        
        int i = 1;
        while (nodes.size() != 0) {
            Node node = nodes.pop();
            node.count = i;
            heap.push(node);
            i++;
        }
        
        return;
    }
    
    for (Node& node: dict.values()) {
        if (node.count > 0) heap.push(node);
    }
}


void Compressor::buildTree()
{
    while (heap.size() > 1) {
        Node left = heap.pop();
        Node right = heap.pop();
        
        if (DefaultLessComparator(right, left)) {
            std::swap(left, right);
        }
        
        left.hasCode = true; left.code = 0;
        right.hasCode = true; right.code = 1;
        Node node;
        
        node.left = new Node(left);
        node.right = new Node(right);
        node.count = left.count + right.count;
        node.hasCode = false;
        
        heap.push(node);
    }
    
    tree = std::shared_ptr<Node>(new Node(heap.pop()));
}


void Compressor::buildMap(std::vector<byte> path, Node* node)
{
    if (node == NULL) node = tree.get();
    
    if (node->hasCode) path.push_back(node->code);
    
    if (node->left != NULL) buildMap(path, node->left);
    map.set(node->key, Bytes(path));
    if (node->right != NULL) buildMap(path, node->right);
    
    if (node->hasCode) path.pop_back();
}


void Compressor::buildHeader(Node* node)
{
    if (node == NULL) node = tree.get();
    
    if (node->left != NULL) {
        header.append("D");
        buildHeader(node->left);
        header.append("U");
    }
    
    if ((node->left == NULL) && (node->right == NULL)) {
        header.append(" " + std::to_string(node->key.length()));
        header.append(" " + node->key + " ");
    }
    
    if (node->right != NULL) {
        buildHeader(node->right);
        header.append("U");
    }
}


void Compressor::writeBits(CBitOStream &os, const Bytes &bytes)
{
    for (const byte& b: bytes.data) os.WriteBit(b);
}


void Compressor::writeHeader(CBitOStream &os)
{
    if (USE_TABLE_HEADER) {
        Heap<Node> nodes;
        
        for (Node& node: dict.values()) {
            if (node.count > 0) nodes.push(node);
        }
        
        while (nodes.size() != 0) {
            Node node = nodes.pop();
            writeByteLen(os, (int) node.key.length());
            
            for (auto& ch: node.key) {
                os.WriteByte(ch);
            }
        }
        
        writeByteLen(os, 0);
        return;
    }
    
    
    int cmdLength = 0;
    Bytes buff;
    
    for (int i = 0; i < (int) header.length(); i++) {
        bool cmdEnd = false;
        
        if ((header[i] == 'U') || (header[i] == 'D')) {
            cmdLength++;
            cmdEnd = true;
            
            switch (header[i]) {
                case 'U': buff.data.push_back(0); break;
                case 'D': buff.data.push_back(1); break;
                default: break;
            }
        }
        
        if ((!cmdEnd) || (i == (int) header.length() - 1)) {
            
            writeLen(os, cmdLength);
            writeBits(os, buff);
            
            if (i == (int) header.length() - 1) {
                writeByteLen(os, 0);
                break;
            }
            
            int len = 0;
            int j = i + 1;
            
            for (; header[j] != ' '; j++) {
                len = len * 10 + (header[j] - 48);
            }
            
            i = j;
            j = i + 1;
            writeByteLen(os, len);
            
            for (int k = 0; k < len; k++, j++) {
                os.WriteByte(header[j]);
            }
            
            i = j;
            cmdLength = 0;
            buff.data.clear();
        }
    }
}


void Compressor::writeLen(CBitOStream &os, const int len)
{
    for (int i = 0; i < len; i++) {
        os.WriteBit(true);
    }
    
    os.WriteBit(false);
}


void Compressor::writeByteLenExplicit(CBitOStream &os, const int len, int bits)
{
    for (int i = 0; i < bits; i++) {
        byte mask = 1 << (bits - 1 - i);
        bool value = (len & mask) != 0;
        os.WriteBit(value);
    }
}


void Compressor::writeByteLen(CBitOStream &os, const int len)
{
    if (!USE_STATIC_LEN) {
        writeLen(os, len);
        return;
    }
    
    writeByteLenExplicit(os, len, LEN_BITS);
}


int Compressor::readByteLenExplicit(CBitIStream &is, int bits)
{
    byte value = 0;
    
    for (int i = 0; i < bits; i++) {
        value = value << 1;
        if (is.ReadBit()) value++;
    }
    
    return value;
}


int Compressor::readByteLen(CBitIStream &is)
{
    if (!USE_STATIC_LEN) return readLen(is);
    return readByteLenExplicit(is, LEN_BITS);
}


Bytes Compressor::readBits(CBitIStream &is, const int len)
{
    Bytes bytes;
    
    for (int i = 0; i < len; i++) {
        bytes.data.push_back((is.ReadBit()) ? (1) : (0));
    }
    
    return bytes;
}


void Compressor::readHeader(CBitIStream &is)
{
    if (USE_TABLE_HEADER) {
        dict.clear();
        int i = 1;
        
        while (!is.IsFinished()) {
            int len = readByteLen(is);
            if (len == 0) break;
            
            std::string key;
            
            for (int i = 0; i < len; i++) {
                key.push_back(is.ReadByte());
            }
            
            Node node = Node(key);
            node.count = i;
            dict.put(key, node);
            i++;
            
            header.append(std::to_string(len) + " ");
            header.append("" + key + " ");
        }
        
        return;
    }
    
    header.clear();
    
    while (!is.IsFinished()) {
        int len = 0;
        len = readLen(is);
        Bytes commands = readBits(is, len);
        
        for (int i = 0; i < len; i++) {
            if (commands.data[i] == 1) {
                header.append("D");
            } else {
                header.append("U");
            }
        }
        
        header.append(" ");
        len = readByteLen(is);
        
        if (len == 0) {
            break;
        }
        
        header.append(std::to_string(len));
        header.append(" ");
        
        for (int i = 0; i < len; i++) {
            header.push_back(is.ReadByte());
        }
        
        header.append(" ");
    }
}


int Compressor::readLen(CBitIStream &is)
{
    int i = 0;
    while (is.ReadBit()) i++;
    return i;
}



void Compressor::parseHeader()
{
    if (USE_TABLE_HEADER) {
        heap.clear();
        buildHeap();
        buildTree();
        return;
    }
    
    Node* current = new Node();
    Node* root = current;
    
    for (int i = 0; i < (int) header.length(); i++) {
        
        if (header[i] == 'D') {
            Node* newNode = new Node();
            newNode->parent = current;
            current->left = newNode;
            newNode->code = 0;
            newNode->hasCode = true;
            current = newNode;
        }
        
        if (header[i] == ' ') {
            if (i == (int) header.length() - 1) break;
            
            int len = 0;
            int j = i + 1;
            
            for (; header[j] != ' '; j++) {
                len = len * 10 + (header[j] - 48);
            }
            
            i = j;
            j = i + 1;
            std::string buff;
            
            for (int k = 0; k < len; k++, j++) {
                buff.push_back(header[j]);
            }
            
            i = j;
            current->key = buff;
        }
        
        if (header[i] == 'U') {
            current = current->parent;
            
            if (current->right == NULL) {
                Node* newNode = new Node();
                newNode->parent = current;
                current->right = newNode;
                newNode->code = 1;
                newNode->hasCode = true;
                current = newNode;
            }
        }
    }
    
    tree = std::shared_ptr<Node>(root);
}



bool Compressor::bitsCompare(const Bytes &a, const Bytes &b)
{
    if (a.data.size() != b.data.size()) return false;
    
    for (int i = 0; i < a.data.size(); i++) {
        if (a.data[i] != b.data[i]) return false;
    }
    
    return true;
}



void Compressor::writeBitsRLE(CBitOStream &os, const Bytes &bytes)
{
    if (!USE_RLE) {
        writeBits(os, bytes);
        return;
    }
    
    static std::vector<Bytes> chain;
    static bool unique = false;
    bool flushing = false;
    if (bytes.data.size() == 0) flushing = true;
    
    if (!flushing && (chain.size() == 0)) {
        chain.push_back(bytes);
        return;
    }
    
    if (!flushing && (chain.size() == 1)) {
        unique = !bitsCompare(chain[chain.size() - 1], bytes);
        chain.push_back(bytes);
        return;
    }
    
    
    if (!flushing && (!bitsCompare(chain[chain.size() - 1], bytes))) {
        if (unique) {
            chain.push_back(bytes);
        } else {
            flushing = true;
        }
    } else {
        if (unique) {
            flushing = true;
        } else {
            chain.push_back(bytes);
        }
    }
    
    if (chain.size() > RLE_MAX_LEN) {
        chain.pop_back();
        flushing = true;
    }
    
    
    if (flushing) {
        if (USE_RLE_DYNAMIC_LENGTH) {
            os.WriteBit(unique);
            writeLen(os, (int) chain.size());
        } else {
            os.WriteBit(unique);
            writeByteLenExplicit(os, (int) chain.size(), RLE_BITS);
        }
        
        if (unique) {
            for (Bytes& bytes: chain) {
                writeBits(os, bytes);
            }
        } else {
            writeBits(os, chain[0]);
        }
        
        chain.clear();
        if (bytes.data.size() != 0) writeBitsRLE(os, bytes);
    }
}


int Compressor::readBitsRLE(CBitIStream &is)
{
    if (!USE_RLE) return -1;
    bool unique = is.ReadBit();
    int signature;
    
    if (USE_RLE_DYNAMIC_LENGTH) {
        signature = readLen(is);
    } else {
        signature = readByteLenExplicit(is, RLE_BITS);
    }
    
    if (unique) signature *= -1;
    return signature;
}






std::string Compressor::compressed()
{
    if (src.length() < 10) return src;
    
    buildStat();
    buildHeap();
    buildTree();
    buildMap();
    buildHeader();
    
    //map.debug();
    
    CBitOStream os;
    writeHeader(os);
    
    std::string word;
    
    for (int i = 0; i < (int) src.length(); i++) {
        
        if (isLetter(src[i])) {
            word.push_back(src[i]);
        }
        
        if ((!isLetter(src[i])) || (i == (int) src.length() - 1) ||
            (word.length() > MAX_WORD)) {
            
            if (word.length() > 0) {
                if (map.has(word)) {
                    writeBits(os, map.get(word));
                } else {
                    for (auto& ch: word) {
                        std::string one;
                        one.push_back(ch);
                        std::cout << "probing_key " << one << std::endl;
                        
                        writeBits(os, map.get(one));
                    }
                }
                
                word.clear();
            }
            
            if (!isLetter(src[i])) {
                std::string one;
                one.push_back(src[i]);
                
                if (map.has(one)) {
                    writeBits(os, map.get(one));
                }
            }
        }
        
        if ((USE_CYRILLIC_PAIRS) && (i < (int) src.length() - 1)) {
            std::string cyrillic;
            
            if (src[i] < 0) {
                cyrillic.push_back(src[i]);
                cyrillic.push_back(src[i + 1]);
                
                if (map.has(cyrillic)) {
                    writeBits(os, map.get(cyrillic));
                    i++;
                    continue;
                }
            }
        }
    }
    
    os.Close();
    return os.GetResult();
}



std::string Compressor::decompressed()
{
    if (src.length() < 10) return src;
    CBitIStream is(src);
    
    readHeader(is);
    parseHeader();
    buildMap();
    
    //map.debug();
    
    std::string result;
    Node* current = tree.get();
    
    while (!is.IsFinished()) {
        if (is.ReadBit()) {
            current = current->right;
        } else {
            current = current->left;
        }
        
        if ((current->left == NULL) && (current->right == NULL)) {
            result.append(current->key);
            current = tree.get();
        }
    }
    
    return result;
}





void compress_string(const std::string &source, std::string &compressed)
{
    Compressor compressor(source);
    compressed = compressor.compressed();
}


void decompress_string(const std::string &compressed, std::string &result)
{
    Compressor decompressor(compressed);
    result = decompressor.decompressed();
}


void autotest(std::string file)
{
    std::ifstream t(file);
    std::string src((std::istreambuf_iterator<char>(t)), (std::istreambuf_iterator<char>()));
    
    Compressor compressor(src);
    std::string compressed = compressor.compressed();
    
    Compressor decompressor(compressed);
    std::string decompressed = decompressor.decompressed();
    
    std::cout << file << std::endl << "\t";
    std::cout << (src == decompressed) << " ";
    
    //std::cout << compressed << std::endl;
    //std::cout << decompressed << std::endl;
    
    float k = compressed.length() / (float)src.length();
    std::cout << k << std::endl;
}


int main(int argc, const char * argv[]) {
    autotest("/Users/coon/Documents/Kernel/src/Alg3_12/Alg3_12/main.cpp");
    autotest("/Users/coon/Documents/Kernel/src/Alg2_51/Alg2_51/main.cpp");
    autotest("/Users/coon/Documents/Kernel/src/Alg2_34/Alg2_34/main.cpp");
    autotest("/Users/coon/Documents/Kernel/src/Alg2_21/Alg2_21/main.cpp");
    autotest("/Users/coon/Documents/Kernel/src/Alg2_12/Alg2_12/main.cpp");
    autotest("/Users/coon/Desktop/test2.txt");
 
 return 0;
}

