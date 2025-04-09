#ifndef LYJ_SWITCHML_HEADER_H
#define LYJ_SWITCHML_HEADER_H

#include "ns3/header.h"

namespace ns3{
/**
 * \ingroup network
 * A simple example of an Header implementation
 */
class SwitchHeader : public Header 
{
public:

  SwitchHeader ();
  virtual ~SwitchHeader ();

  /**
   * Set the header data.
   * \param data The data.
   */
  void SetKey (uint32_t data);
  /**
   * Get the header data.
   * \return The data.
   */
  uint32_t GetKey (void) const;

  /**
   * Set the header data.
   * \param data The data.
   */
  void SetHostid (uint8_t data);
  /**
   * Get the header data.
   * \return The data.
   */
  uint8_t GetHostid (void) const;

  void SetHostnum (uint8_t data);
  /**
   * Get the header data.
   * \return The data.
   */
  uint8_t GetHostnum (void) const;

  void SetDelta (uint32_t tmpdelta);
  /**
   * Get the header data.
   * \return The data.
   */
  uint32_t GetDelta (void) const;

   /**
   * Set the header data.
   * \param data The data.
   */
  void SetACK (uint8_t data);
  /**
   * Get the header data.
   * \return The data.
   */
  uint8_t GetACK (void) const;

  void SetMACK (uint8_t data);
  uint8_t GetMACK (void) const;


  void SetMerged(uint8_t isMerged);
  uint8_t GetMerged (void) const;

  void SetTotal(uint16_t m_total);
  uint16_t GetTotal (void) const;

  void SetbpAggr (uint8_t data);
  uint8_t GetbpAggr (void) const;

  void SetCollision (uint16_t data);
  /**
   * Get the header data.
   * \return The data.
   */
  uint16_t GetCollision (void) const;

  // void SetEnd (uint16_t data);
  // /**
  //  * Get the header data.
  //  * \return The data.
  //  */
  // uint16_t GetEnd (void) const;

  void SetAppID (uint16_t data);
  /**
   * Get the header data.
   * \return The data.
   */
  uint16_t GetAppID (void) const;


  /**
   * \brief Get the type ID.
   * \return the object TypeId
   */
  static TypeId GetTypeId (void);
  virtual TypeId GetInstanceTypeId (void) const;
  virtual void Print (std::ostream &os) const;
  virtual void Serialize (Buffer::Iterator start) const;
  virtual uint32_t Deserialize (Buffer::Iterator start);
  virtual uint32_t GetSerializedSize (void) const;
private:
  uint32_t key;  //!< Header data
  uint8_t host_id;
  uint8_t host_num; 
  uint8_t isAck; //
  uint8_t isMACK; //lzy
  uint8_t bpAggr; //for a2tp
  uint16_t isCollision;
  uint16_t m_total;
  uint8_t isMerged; // wyq
  //uint16_t isEnd;
  uint16_t app_id;
  uint32_t delta;
};
}

#endif /* NS3_SOCKET_H */