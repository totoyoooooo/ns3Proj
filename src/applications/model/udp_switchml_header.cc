#include "udp_switchml_header.h"

namespace ns3{
    
SwitchHeader::SwitchHeader ()
{
  // we must provide a public default constructor, 
  // implicit or explicit, but never private.
}
SwitchHeader::~SwitchHeader ()
{
}

TypeId
SwitchHeader::GetTypeId (void)
{
  static TypeId tid = TypeId ("ns3::SwitchHeader")
    .SetParent<Header> ()
    .AddConstructor<SwitchHeader> ()
  ;
  return tid;
}
TypeId
SwitchHeader::GetInstanceTypeId (void) const
{
  return GetTypeId ();
}

void
SwitchHeader::Print (std::ostream &os) const
{
  // This method is invoked by the packet printing
  // routines to print the content of my header.
  //os << "data=" << m_data << std::endl;
  os << "data=" << key;
}

void 
SwitchHeader::SetMerged (uint8_t ismerged)
{
  isMerged = ismerged;
}
uint8_t 
SwitchHeader::GetMerged (void) const
{
  return isMerged;
}
void 
SwitchHeader::SetTotal (uint16_t total)
{
  m_total = total;
}
uint16_t 
SwitchHeader::GetTotal (void) const
{
  return m_total;
}
uint32_t
SwitchHeader::GetSerializedSize (void) const
{
  // we reserve 2 bytes for our header.
  return 20; // wyq （新增3字节）
}
void
SwitchHeader::Serialize (Buffer::Iterator start) const
{
  // we can serialize two bytes at the start of the buffer.
  // we write them in network byte order.
  start.WriteHtonU32 (key);
  //start.WriteHtonU16 (host_id);
  start.WriteU8(host_id);
  start.WriteU8(host_num);
  start.WriteU8(isAck);
  start.WriteU8(bpAggr);
  start.WriteU8(isMACK);
  start.WriteHtonU16(m_total);
  start.WriteU8(isMerged);
  start.WriteHtonU16 (isCollision);
  //start.WriteHtonU16 (isEnd);
  start.WriteHtonU16(app_id);
  start.WriteHtonU32 (delta);
}
uint32_t
SwitchHeader::Deserialize (Buffer::Iterator start)
{
  // we can deserialize two bytes from the start of the buffer.
  // we read them in network byte order and store them
  // in host byte order.
  key = start.ReadNtohU32 ();
  //host_id = start.ReadNtohU16 ();
  host_id = start.ReadU8();
  host_num = start.ReadU8();
  //isAck = start.ReadNtohU16();
  isAck = start.ReadU8();
  bpAggr = start.ReadU8();
  isMACK = start.ReadU8();
  m_total = start.ReadNtohU16();
  isMerged = start.ReadU8();
  isCollision = start.ReadNtohU16();
  //isEnd = start.ReadNtohU16();
  app_id = start.ReadNtohU16();
  delta = start.ReadNtohU32();
  // we return the number of bytes effectively read.
  return 20;
}

void 
SwitchHeader::SetKey (uint32_t data)
{
  key = data;
}
uint32_t 
SwitchHeader::GetKey (void) const
{
  return key;
}

//lzy
void 
SwitchHeader::SetDelta (uint32_t tmpdelta)
{
  delta = tmpdelta;
}
uint32_t 
SwitchHeader::GetDelta (void) const
{
  return delta;
}

void 
SwitchHeader::SetHostid (uint8_t data)
{
  host_id = data;
}
uint8_t 
SwitchHeader::GetHostid (void) const
{
  return host_id;
}

void 
SwitchHeader::SetHostnum (uint8_t data)
{
  host_num = data;
}
uint8_t 
SwitchHeader::GetHostnum (void) const
{
  return host_num;
}

void 
SwitchHeader::SetACK (uint8_t data)
{
  isAck = data;
}
uint8_t 
SwitchHeader::GetACK (void) const
{
  return isAck;
}

void 
SwitchHeader::SetMACK (uint8_t data)
{
  isMACK = data;
}
uint8_t 
SwitchHeader::GetMACK (void) const
{
  return isMACK;
}

void 
SwitchHeader::SetbpAggr (uint8_t data)
{
  bpAggr = data;
}
uint8_t 
SwitchHeader::GetbpAggr (void) const
{
  return bpAggr;
}

void 
SwitchHeader::SetCollision (uint16_t data)
{
  isCollision = data;
}
uint16_t 
SwitchHeader::GetCollision (void) const
{
  return isCollision;
}

// void 
// SwitchHeader::SetEnd (uint16_t data)
// {
//   isEnd = data;
// }
// uint16_t 
// SwitchHeader::GetEnd (void) const
// {
//   return isEnd;
// }

void 
SwitchHeader::SetAppID (uint16_t data)
{
  app_id = data;
}
uint16_t 
SwitchHeader::GetAppID (void) const
{
  return app_id;
}

}