console.log("here is the call");
var url = 'http://127.0.0.1:8000/purchase/page1/result/orders/';
var data =
{
    "id": "4",
    "order_first_name": "Robert",
    "order_last_name": "Ohr",
    "order_email": "bob@tca.io",
    "order_wallet_address": "bobrichwallet",
    "order_street_address": "Boston Cambridge",
    "order_price": "8800",
    "order_completed": "incomplete",
    "pub_date": "2021-09-13T23:02:39Z"
};


async function postData()
{
  const response = await fetch(url,
  {
    method: 'POST',
    mode:'cors',
    headers:
    {
      'Content-type': 'application/json',
    },
    body: JSON.stringify(data),
  })
  .then(response => response.json())
  .then(data =>
    {
      console.log('Success:', data);
      console.log("response: "+response);
    })
  .catch((error) => { console.error('Error:', error); });
}
postData();
/*
var request = new XMLHttpRequest();
request.open('POST', 'http://127.0.0.1:8000/purchase/page1/result/orders/', true);
request.setRequestHeader("Content-type", "application/json");
request.send(data);
*/
