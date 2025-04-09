/* -*- Mode:C++; c-file-style:"gnu"; indent-tabs-mode:nil; -*- */
/*
 * Copyright (c) 2008 INRIA
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License version 2 as
 * published by the Free Software Foundation;
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 * Author: Mathieu Lacage <mathieu.lacage@sophia.inria.fr>
 */
#include "udp-worker-helper.h"
#include "ns3/udp-worker.h"
#include "ns3/udp-aggregator.h"
#include "ns3/udp-ps.h"
#include "ns3/uinteger.h"
#include "ns3/names.h"

namespace ns3 {

UdpPsHelper::UdpPsHelper (uint16_t port)
{
  m_factory.SetTypeId (UdpPs::GetTypeId ());
  SetAttribute ("Port", UintegerValue (port));
}

void 
UdpPsHelper::SetAttribute (
  std::string name, 
  const AttributeValue &value)
{
  m_factory.Set (name, value);
}

ApplicationContainer
UdpPsHelper::Install (Ptr<Node> node) const
{
  return ApplicationContainer (InstallPriv (node));
}

ApplicationContainer
UdpPsHelper::Install (std::string nodeName) const
{
  Ptr<Node> node = Names::Find<Node> (nodeName);
  return ApplicationContainer (InstallPriv (node));
}

ApplicationContainer
UdpPsHelper::Install (NodeContainer c) const
{
  ApplicationContainer apps;
  for (NodeContainer::Iterator i = c.Begin (); i != c.End (); ++i)
    {
      apps.Add (InstallPriv (*i));
    }

  return apps;
}

Ptr<Application>
UdpPsHelper::InstallPriv (Ptr<Node> node) const
{
  Ptr<Application> app = m_factory.Create<UdpPs> ();
  node->AddApplication (app);

  return app;
}

UdpAggregatorHelper::UdpAggregatorHelper (uint16_t port)
{
  m_factory.SetTypeId (UdpAggregator::GetTypeId ());
  SetAttribute ("Port", UintegerValue (port));
}

void 
UdpAggregatorHelper::SetAttribute (
  std::string name, 
  const AttributeValue &value)
{
  m_factory.Set (name, value);
}

ApplicationContainer
UdpAggregatorHelper::Install (Ptr<Node> node) const
{
  return ApplicationContainer (InstallPriv (node));
}

ApplicationContainer
UdpAggregatorHelper::Install (std::string nodeName) const
{
  Ptr<Node> node = Names::Find<Node> (nodeName);
  return ApplicationContainer (InstallPriv (node));
}

ApplicationContainer
UdpAggregatorHelper::Install (NodeContainer c) const
{
  ApplicationContainer apps;
  for (NodeContainer::Iterator i = c.Begin (); i != c.End (); ++i)
    {
      apps.Add (InstallPriv (*i));
    }

  return apps;
}

Ptr<Application>
UdpAggregatorHelper::InstallPriv (Ptr<Node> node) const
{
  Ptr<Application> app = m_factory.Create<UdpAggregator> ();
  node->AddApplication (app);

  return app;
}


UdpWorkerHelper::UdpWorkerHelper (Address address, uint16_t port)
{
  m_factory.SetTypeId (UdpWorker::GetTypeId ());
  SetAttribute ("RemoteAddress", AddressValue (address));
  SetAttribute ("RemotePort", UintegerValue (port));
}

UdpWorkerHelper::UdpWorkerHelper (Address address)
{
  m_factory.SetTypeId (UdpWorker::GetTypeId ());
  SetAttribute ("RemoteAddress", AddressValue (address));
}

void 
UdpWorkerHelper::SetAttribute (
  std::string name, 
  const AttributeValue &value)
{
  m_factory.Set (name, value);
}

void
UdpWorkerHelper::SetFill (Ptr<Application> app, std::string fill)
{
  app->GetObject<UdpWorker>()->SetFill (fill);
}

void
UdpWorkerHelper::SetFill (Ptr<Application> app, uint8_t fill, uint32_t dataLength)
{
  app->GetObject<UdpWorker>()->SetFill (fill, dataLength);
}

void
UdpWorkerHelper::SetFill (Ptr<Application> app, uint8_t *fill, uint32_t fillLength, uint32_t dataLength)
{
  app->GetObject<UdpWorker>()->SetFill (fill, fillLength, dataLength);
}

ApplicationContainer
UdpWorkerHelper::Install (Ptr<Node> node) const
{
  return ApplicationContainer (InstallPriv (node));
}

ApplicationContainer
UdpWorkerHelper::Install (std::string nodeName) const
{
  Ptr<Node> node = Names::Find<Node> (nodeName);
  return ApplicationContainer (InstallPriv (node));
}

ApplicationContainer
UdpWorkerHelper::Install (NodeContainer c) const
{
  ApplicationContainer apps;
  for (NodeContainer::Iterator i = c.Begin (); i != c.End (); ++i)
    {
      apps.Add (InstallPriv (*i));
    }

  return apps;
}

Ptr<Application>
UdpWorkerHelper::InstallPriv (Ptr<Node> node) const
{
  Ptr<Application> app = m_factory.Create<UdpWorker> ();
  node->AddApplication (app);

  return app;
}

} // namespace ns3
